#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfFormat.h>

#include <pbbam/StringUtilities.h>
#include <pbbam/vcf/VcfHeader.h>
#include "VcfFormatException.h"

#include <htslib/vcf.h>

#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <cassert>
#include <cmath>

namespace PacBio {
namespace VCF {
namespace {

// using htslib's current version for better compatibility
static constexpr const char current_version[] = "VCFv4.2";

namespace Tokens {

static constexpr const char file_format[] = "fileformat";

static constexpr const char double_hash[] = "##";
static constexpr const char contig_lead[] = "##contig=<";
static constexpr const char filter_lead[] = "##FILTER=<";
static constexpr const char format_lead[] = "##FORMAT=<";
static constexpr const char info_lead[] = "##INFO=<";
static constexpr const char chrom_lead[] = "#CHROM";

static constexpr const char id[] = "ID";
static constexpr const char number[] = "Number";
static constexpr const char type[] = "Type";
static constexpr const char description[] = "Description";
static constexpr const char source[] = "Source";
static constexpr const char version[] = "Version";

}  // namespace Tokens

std::string QuotedText(const std::string& d) { return "\"" + d + "\""; }

std::string UnquotedText(const std::string& d)
{
    if (d.size() < 2 || d.front() != '"' || d.back() != '"') {
        throw VcfFormatException{"description text is not quoted: " + d};
    }
    return d.substr(1, d.size() - 2);
}

}  // namespace

const char* VcfFormat::CurrentVersion() { return current_version; }

std::string VcfFormat::FormattedContigDefinition(const ContigDefinition& def)
{
    std::ostringstream text;

    // ID
    text << Tokens::contig_lead << Tokens::id << '=' << def.Id();

    // attributes
    if (!def.Attributes().empty()) {
        text << ',';
        bool first = true;
        for (const auto& attr : def.Attributes()) {
            if (!first) {
                text << ',';
            }
            text << attr.first << '=' << attr.second;
            first = false;
        }
    }
    text << '>';
    return text.str();
}

std::string VcfFormat::FormattedFilterDefinition(const FilterDefinition& def)
{
    std::ostringstream text;
    text << Tokens::filter_lead << Tokens::id << '=' << def.Id() << ',' << Tokens::description
         << '=' << QuotedText(def.Description()) << '>';
    return text.str();
}

std::string VcfFormat::FormattedFormatDefinition(const FormatDefinition& def)
{
    std::ostringstream text;
    text << Tokens::format_lead << Tokens::id << '=' << def.Id() << ',' << Tokens::number << '='
         << def.Number() << ',' << Tokens::type << '=' << def.Type() << ',' << Tokens::description
         << '=' << QuotedText(def.Description()) << '>';
    return text.str();
}

std::string VcfFormat::FormattedGeneralDefinition(const GeneralDefinition& def)
{
    std::ostringstream text;
    text << Tokens::double_hash << def.Id() << '=' << def.Text();
    return text.str();
}

std::string VcfFormat::FormattedInfoDefinition(const InfoDefinition& def)
{
    std::ostringstream text;
    text << Tokens::info_lead << Tokens::id << '=' << def.Id() << ',' << Tokens::number << '='
         << def.Number() << ',' << Tokens::type << '=' << def.Type() << ',' << Tokens::description
         << '=' << QuotedText(def.Description());

    if (def.Source() && !def.Source()->empty()) {
        text << ',' << Tokens::source << '=' << QuotedText(*def.Source());
    }

    if (def.Version() && !def.Version()->empty()) {
        text << ',' << Tokens::version << '=' << QuotedText(*def.Version());
    }

    text << '>';
    return text.str();
}

std::string VcfFormat::FormattedHeader(const VcfHeader& header)
{
    std::ostringstream out;

    const auto& fileformat = header.GeneralDefinition(Tokens::file_format);
    out << FormattedGeneralDefinition(fileformat) << '\n';

    // remaining general definiitions
    for (const auto& def : header.GeneralDefinitions()) {
        if (def.Id() != Tokens::file_format) {
            out << FormattedGeneralDefinition(def) << '\n';
        }
    }

    // ##contig
    for (const auto& contig : header.ContigDefinitions()) {
        out << FormattedContigDefinition(contig) << '\n';
    }

    // ##FILTER
    for (const auto& filter : header.FilterDefinitions()) {
        out << FormattedFilterDefinition(filter) << '\n';
    }

    // ##INFO
    for (const auto& info : header.InfoDefinitions()) {
        out << FormattedInfoDefinition(info) << '\n';
    }

    // ##FORMAT
    for (const auto& format : header.FormatDefinitions()) {
        out << FormattedFormatDefinition(format) << '\n';
    }

    // fixed headers
    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

    // samples
    const auto& samples = header.Samples();
    if (!samples.empty()) {
        out << "\tFORMAT";
        for (const auto& sample : samples) {
            out << '\t' << sample;
        }
    }

    return out.str();
}

ContigDefinition VcfFormat::ParsedContigDefinition(std::string line)
{
    // should already be checked by "normal" code path
    assert(line.find(Tokens::contig_lead) == 0);

    // substring between brackets
    const auto lastBracketPos = line.find_last_of('>');
    if (lastBracketPos == std::string::npos) {
        throw VcfFormatException{"malformed ##contig line: " + line};
    }
    line = std::string(line.cbegin() + 10, line.cbegin() + lastBracketPos);

    std::string id;
    std::vector<std::pair<std::string, std::string>> attributes;

    const auto fields = BAM::Split(line, ',');
    for (const auto& field : fields) {
        const auto tokens = BAM::Split(field, '=');
        if (tokens.size() != 2) {
            throw VcfFormatException{"malformed ##contig line: " + line};
        }

        if (tokens[0] == Tokens::id) {
            id = tokens[1];
        } else {
            attributes.push_back(std::make_pair(tokens[0], tokens[1]));
        }
    }

    return ContigDefinition{std::move(id), std::move(attributes)};
}

FilterDefinition VcfFormat::ParsedFilterDefinition(std::string line)
{
    // should already be checked by "normal" code path
    assert(line.find(Tokens::filter_lead) == 0);

    // substring between brackets
    const auto lastBracketPos = line.find_last_of('>');
    if (lastBracketPos == std::string::npos) {
        throw VcfFormatException{"malformed FILTER line: " + line};
    }
    line = std::string(line.cbegin() + 10, line.cbegin() + lastBracketPos);

    std::string id;
    std::string description;

    const auto fields = BAM::Split(line, ',');
    for (const auto& field : fields) {
        const auto tokens = BAM::Split(field, '=');
        if (tokens.size() != 2) {
            throw VcfFormatException{"malformed FILTER line: " + line};
        }

        if (tokens[0] == Tokens::id) {
            id = tokens[1];
        } else if (tokens[0] == Tokens::description) {
            description = UnquotedText(tokens[1]);
        } else {
            throw VcfFormatException{"unrecognized FILTER field: " + tokens[0]};
        }
    }

    return FilterDefinition{std::move(id), std::move(description)};
}

FormatDefinition VcfFormat::ParsedFormatDefinition(std::string line)
{
    // should already be checked by "normal" code path
    assert(line.find(Tokens::format_lead) == 0);

    // substring between brackets
    const auto lastBracketPos = line.find_last_of('>');
    if (lastBracketPos == std::string::npos) {
        throw VcfFormatException{"malformed FORMAT line: " + line};
    }
    line = std::string(line.cbegin() + 10, line.cbegin() + lastBracketPos);

    std::string id;
    std::string number;
    std::string type;
    std::string description;

    const auto fields = BAM::Split(line, ',');
    for (const auto& field : fields) {
        const auto tokens = BAM::Split(field, '=');
        if (tokens.size() != 2) {
            throw VcfFormatException{"malformed FORMAT line: " + line};
        }
        if (tokens[0] == Tokens::id) {
            id = tokens[1];
        } else if (tokens[0] == Tokens::number) {
            number = tokens[1];
        } else if (tokens[0] == Tokens::type) {
            type = tokens[1];
        } else if (tokens[0] == Tokens::description) {
            description = UnquotedText(tokens[1]);
        } else {
            throw VcfFormatException{"unrecognized FORMAT field: " + tokens[0]};
        }
    }

    return FormatDefinition{std::move(id), std::move(number), std::move(type),
                            std::move(description)};
}

GeneralDefinition VcfFormat::ParsedGeneralDefinition(const std::string& line)
{
    const auto tokens = BAM::Split(line, '=');
    if (tokens.size() != 2 || tokens[0].find(Tokens::double_hash) != 0) {
        throw VcfFormatException{"malformed header line: " + line};
    }
    return GeneralDefinition{tokens[0].substr(2), tokens[1]};
}

InfoDefinition VcfFormat::ParsedInfoDefinition(std::string line)
{
    // should already be checked by "normal" code path
    assert(line.find(Tokens::info_lead) == 0);

    // substring between brackets
    const auto lastBracketPos = line.find_last_of('>');
    if (lastBracketPos == std::string::npos) {
        throw VcfFormatException{"malformed INFO line: " + line};
    }
    line = std::string(line.cbegin() + 8, line.cbegin() + lastBracketPos);

    std::string id;
    std::string number;
    std::string type;
    std::string description;
    std::string source;
    std::string version;

    const auto fields = BAM::Split(line, ',');
    for (const auto& field : fields) {
        const auto tokens = BAM::Split(field, '=');
        if (tokens.size() != 2) {
            throw VcfFormatException{"malformed INFO line: " + line};
        }

        if (tokens[0] == Tokens::id) {
            id = tokens[1];
        } else if (tokens[0] == Tokens::number) {
            number = tokens[1];
        } else if (tokens[0] == Tokens::type) {
            type = tokens[1];
        } else if (tokens[0] == Tokens::description) {
            description = UnquotedText(tokens[1]);
        } else if (tokens[0] == Tokens::source) {
            source = UnquotedText(tokens[1]);
        } else if (tokens[0] == Tokens::version) {
            version = UnquotedText(tokens[1]);
        } else {
            throw VcfFormatException{"unrecognized INFO field: " + tokens[0]};
        }
    }

    return InfoDefinition{std::move(id),          std::move(number), std::move(type),
                          std::move(description), std::move(source), std::move(version)};
}

VcfHeader VcfFormat::ParsedHeader(const std::string& hdrText)
{
    VcfHeader hdr;

    std::istringstream text{hdrText};
    std::string line;

    // quick check for fileformat - should be the first line
    std::getline(text, line);
    auto genDef = ParsedGeneralDefinition(line);
    if (genDef.Id() != Tokens::file_format) {
        throw VcfFormatException{"file must begin with #fileformat line"};
    }
    hdr.AddGeneralDefinition(std::move(genDef));

    // read through rest of header
    bool chromLineFound = false;
    for (; std::getline(text, line);) {
        if (line.empty()) {
            continue;
        }

        // info line
        if (line.find(Tokens::info_lead) == 0) {
            hdr.AddInfoDefinition(ParsedInfoDefinition(line));

            // filter line
        } else if (line.find(Tokens::filter_lead) == 0) {
            hdr.AddFilterDefinition(ParsedFilterDefinition(line));

            // format line
        } else if (line.find(Tokens::format_lead) == 0) {
            hdr.AddFormatDefinition(ParsedFormatDefinition(line));

            // contig line
        } else if (line.find(Tokens::contig_lead) == 0) {
            hdr.AddContigDefinition(ParsedContigDefinition(line));

            // general comment line
            //
            // NOTE: Check this after all other specific header line types. This
            //       catches all remaining lines starting with "##"
            //
        } else if (line.find(Tokens::double_hash) == 0) {
            hdr.AddGeneralDefinition(ParsedGeneralDefinition(line));

            // CHROM line (maybe w/ samples)
        } else if (line.find(Tokens::chrom_lead) == 0) {
            std::vector<Sample> samples;

            // If samples are present, skip the fixed colums & FORMAT column (9)
            // and read the remaining column labels as sample names.
            //
            auto columns = BAM::Split(line, '\t');
            for (size_t i = 9; i < columns.size(); ++i) {
                samples.push_back(std::move(columns[i]));
            }
            hdr.Samples(std::move(samples));

            // quit header parsing after CHROM line
            chromLineFound = true;
            break;
        } else {
            throw VcfFormatException{"unexpected line found in header:\n" + line};
        }
    }

    if (!chromLineFound) {
        throw VcfFormatException{"CHROM column line is missing"};
    }
    return hdr;
}

VcfHeader VcfFormat::HeaderFromFile(const std::string& fn)
{
    std::ifstream in(fn);
    return HeaderFromStream(in);
}

VcfHeader VcfFormat::HeaderFromStream(std::istream& in)
{
    std::ostringstream text;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line.front() == '#') {
            text << line << '\n';
        } else {
            break;
        }
    }

    return ParsedHeader(text.str());
}

InfoField VcfFormat::ParsedInfoField(const std::string& text)
{
    const auto tokens = BAM::Split(text, '=');
    if (tokens.empty()) {
        throw VcfFormatException{"malformed INFO field: " + text};
    }

    // required ID
    InfoField result;
    result.id = tokens.at(0);
    if (tokens.size() == 1) {
        return result;
    }

    // optional value or values
    const auto& valueStr = tokens.at(1);
    const auto commaFound = valueStr.find(',');
    if (commaFound != std::string::npos) {
        result.values = BAM::Split(valueStr, ',');
    } else {
        result.value = valueStr;
    }

    return result;
}

std::vector<InfoField> VcfFormat::ParsedInfoFields(const std::string& text)
{
    std::vector<InfoField> result;
    const auto& fields = BAM::Split(text, ';');
    for (const auto& field : fields) {
        result.push_back(ParsedInfoField(field));
    }
    return result;
}

GenotypeField VcfFormat::ParsedGenotypeField(const std::string& field)
{
    GenotypeField result;
    const auto fieldValues = BAM::Split(field, ':');
    for (const auto& fieldValue : fieldValues) {
        GenotypeData data;
        auto genotypeDataValues = BAM::Split(fieldValue, ',');
        if (genotypeDataValues.size() == 1) {
            data.value = std::move(genotypeDataValues[0]);
        } else {
            data.values = std::move(genotypeDataValues);
        }
        result.data.push_back(std::move(data));
    }
    return result;
}

VcfVariant VcfFormat::ParsedVariant(const std::string& line)
{
    auto fields = BAM::Split(line, '\t');
    if (fields.size() < 7) {
        throw VcfFormatException{"record is missing required fields: " + line};
    }

    // CHROM POS ID REF ALT REF
    auto chrom = std::move(fields[0]);
    auto pos = std::stoi(fields[1]);
    auto id = std::move(fields[2]);
    auto ref = std::move(fields[3]);
    auto alt = std::move(fields[4]);

    VcfVariant var{std::move(id), std::move(chrom), std::move(pos), std::move(ref), std::move(alt)};

    // QUAL
    const auto& qualStr = fields[5];
    const float qual = (qualStr == "." ? NAN : std::stof(qualStr));
    var.Quality(qual);

    // FILTER
    var.Filter(std::move(fields[6]));

    // INFO (allow empty)
    if (fields.size() >= 8) {
        var.InfoFields(ParsedInfoFields(fields.at(7)));
    }

    // GENOTYPE (samples)
    if (fields.size() > 9) {
        var.GenotypeIds(BAM::Split(fields.at(8), ':'));

        std::vector<GenotypeField> genotypes;
        for (size_t i = 9; i < fields.size(); ++i) {
            genotypes.emplace_back(ParsedGenotypeField(fields.at(i)));
        }
        var.Genotypes(std::move(genotypes));
    }

    return var;
}

std::string VcfFormat::FormattedInfoField(const InfoField& field)
{
    std::ostringstream out;
    out << field.id;
    if (field.value) {
        out << '=' << *field.value;
    } else if (field.values) {
        out << '=' << BAM::Join(*field.values, ',');
    }
    return out.str();
}

std::string VcfFormat::FormattedInfoFields(const std::vector<InfoField>& fields)
{
    std::vector<std::string> result;
    for (const auto& field : fields) {
        result.push_back(FormattedInfoField(field));
    }
    return BAM::Join(result, ';');
}

std::string VcfFormat::FormattedGenotypeField(const GenotypeField& field)
{
    std::string result;
    bool firstDataEntry = true;
    for (const auto& d : field.data) {
        if (!firstDataEntry) {
            result += ':';
        }
        if (d.value) {
            result += *d.value;
        } else {
            assert(d.values);
            result += BAM::Join(*d.values, ',');
        }
        firstDataEntry = false;
    }
    return result;
}

std::string VcfFormat::FormattedVariant(const VcfVariant& var)
{
    std::ostringstream out;
    out << var.Chrom() << '\t' << var.Position() << '\t' << var.Id() << '\t' << var.RefAllele()
        << '\t' << var.AltAllele() << '\t'
        << (var.IsQualityMissing() ? "." : std::to_string(var.Quality())) << '\t' << var.Filter()
        << '\t' << FormattedInfoFields(var.InfoFields());

    const auto genotypeIds = var.GenotypeIds();
    if (!genotypeIds.empty()) {
        out << '\t' << BAM::Join(genotypeIds, ':');
        const auto genotypes = var.Genotypes();
        for (const auto& genotype : genotypes) {
            out << '\t' << FormattedGenotypeField(genotype);
        }
    }
    return out.str();
}

}  // namespace VCF
}  // namespace PacBio

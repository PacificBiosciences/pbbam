// File Description
/// \file CCSRecordFormat.cpp
/// \brief Implements the CCSRecordFormat class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ccs/CCSRecordFormat.h"

#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "pbbam/StringUtilities.h"

namespace {

static const std::string MovieName{"movie_name"};
static const std::string BindingKit{"binding_kit"};
static const std::string SequencingKit{"sequencing_kit"};
static const std::string BasecallerVersion{"basecaller_version"};
static const std::string FrameRate{"framerate"};

}  // namespace

namespace PacBio {
namespace CCS {

CCSHeader CCSRecordFormat::DeserializeHeader(const std::vector<std::string>& lines)
{
    if (lines.empty())
        throw std::runtime_error{
            "[pbbam] CCS record format ERROR: cannot create header from empty text"};

    CCSHeader result;
    std::vector<std::string> fields;
    for (const auto& line : lines) {
        fields = PacBio::BAM::Split(line, '=');
        if (fields.size() != 2) {
            std::ostringstream msg;
            msg << "[pbbam] CCS record format ERROR: must have syntax 'name=value' in header "
                   "line:\n"
                << line;
            throw std::runtime_error{msg.str()};
        }

        // clang-format off
        if      (fields[0] == MovieName)         result.MovieName = fields[1];
        else if (fields[0] == BindingKit)        result.BindingKit = fields[1];
        else if (fields[0] == SequencingKit)     result.SequencingKit = fields[1];
        else if (fields[0] == BasecallerVersion) result.BasecallerVersion = fields[1];
        else if (fields[0] == FrameRate)         result.FrameRate = fields[1];
        else {
            std::ostringstream msg;
            msg << "[pbbam] CCS record format ERROR: unrecognized header field name: '" << fields[0] << '\'';
            throw std::runtime_error{msg.str()};
        }
        // clang-format on
    }
    return result;
}

std::vector<std::string> CCSRecordFormat::SerializeHeader(const CCSHeader& header)
{
    // clang-format off
    std::vector<std::string> result;
    result.push_back(MovieName         + '=' + header.MovieName);
    result.push_back(BindingKit        + '=' + header.BindingKit);
    result.push_back(SequencingKit     + '=' + header.SequencingKit);
    result.push_back(BasecallerVersion + '=' + header.BasecallerVersion);
    result.push_back(FrameRate         + '=' + header.FrameRate);
    // clang-format on
    return result;
}

CCSRecord CCSRecordFormat::DeserializeRecord(const std::string& line)
{
    const auto fields = PacBio::BAM::Split(line);
    if (fields.size() != 8) {
        std::ostringstream msg;
        msg << "[pbbam] CCS record format ERROR: malformed record line:\n" << line;
        throw std::runtime_error{msg.str()};
    }

    // clang-format off
    CCSRecord result;
    result.HoleNumber = std::stoi(fields[0]);
    result.QueryStart = std::stoi(fields[1]);
    result.QueryEnd   = std::stoi(fields[2]);
    result.LocalContextFlags = static_cast<PacBio::BAM::LocalContextFlags>(std::stoul(fields[3]));
    result.Accuracy = std::stof(fields[4]);

    const auto snrs = PacBio::BAM::Split(fields[5], ',');
    if (snrs.size() != 4) {
        std::ostringstream msg;
        msg << "[pbbam] CCS record format ERROR: SNR field must have 4 values";
        throw std::runtime_error{msg.str()};
    }
    result.SignalToNoise = {
        std::stod(snrs[0]),
        std::stod(snrs[1]),
        std::stod(snrs[2]),
        std::stod(snrs[3])
    };

    result.Sequence = fields[6];

    const auto pwStrings = PacBio::BAM::Split(fields[7], ',');
    std::vector<uint16_t> pws;
    pws.reserve(pwStrings.size());
    for (const auto& pwString : pwStrings)
        pws.emplace_back(std::stoul(pwString));
    result.PulseWidths = pws;

    // clang-format on
    return result;
}

std::string CCSRecordFormat::SerializeRecord(const CCSRecord& record)
{
    // clang-format off
    std::ostringstream out;
    out << record.HoleNumber << '\t'
        << record.QueryStart << '\t'
        << record.QueryEnd << '\t'
        << static_cast<uint16_t>(record.LocalContextFlags) << '\t'
        << record.Accuracy << '\t'
        << record.SignalToNoise.A << ',' << record.SignalToNoise.C << ','
        << record.SignalToNoise.G << ',' << record.SignalToNoise.T << '\t'
        << record.Sequence << '\t';

    bool first = true;
    for (const auto pw : record.PulseWidths) {
        if (!first) out << ',';
        else first = false;
        out << pw;
    }

    // clang-format on
    return out.str();
}

}  // namespace CCS
}  // namespace PacBio
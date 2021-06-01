#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfHeader.h>

#include <type_traits>

#include <pbbam/vcf/VcfFormat.h>

namespace PacBio {
namespace VCF {

VcfHeader::VcfHeader() { Version(VcfFormat::CurrentVersion()); }

VcfHeader::VcfHeader(const std::string& hdrText) { *this = VcfFormat::ParsedHeader(hdrText); }

VcfHeader& VcfHeader::AddContigDefinition(VCF::ContigDefinition contig)
{
    const auto found = contigLookup_.find(contig.Id());
    if (found == contigLookup_.cend()) {
        contigLookup_.insert({contig.Id(), contigDefinitions_.size()});
        contigDefinitions_.push_back(std::move(contig));
    } else {
        contigDefinitions_.at(found->second) = std::move(contig);
    }
    return *this;
}

VcfHeader& VcfHeader::AddFilterDefinition(VCF::FilterDefinition filter)
{
    const auto found = filterLookup_.find(filter.Id());
    if (found == filterLookup_.cend()) {
        filterLookup_.insert({filter.Id(), filterDefinitions_.size()});
        filterDefinitions_.push_back(std::move(filter));
    } else {
        filterDefinitions_.at(found->second) = std::move(filter);
    }
    return *this;
}

VcfHeader& VcfHeader::AddFormatDefinition(VCF::FormatDefinition format)
{
    const auto found = formatLookup_.find(format.Id());
    if (found == formatLookup_.cend()) {
        formatLookup_.insert({format.Id(), formatDefinitions_.size()});
        formatDefinitions_.push_back(std::move(format));
    } else {
        formatDefinitions_.at(found->second) = std::move(format);
    }
    return *this;
}

VcfHeader& VcfHeader::AddGeneralDefinition(VCF::GeneralDefinition def)
{
    const auto found = generalLookup_.find(def.Id());
    if (found == generalLookup_.cend()) {
        generalLookup_.insert({def.Id(), generalDefinitions_.size()});
        generalDefinitions_.push_back(std::move(def));
    } else {
        generalDefinitions_.at(found->second) = std::move(def);
    }
    return *this;
}

VcfHeader& VcfHeader::AddInfoDefinition(VCF::InfoDefinition info)
{
    const auto found = infoLookup_.find(info.Id());
    if (found == infoLookup_.cend()) {
        infoLookup_.insert({info.Id(), infoDefinitions_.size()});
        infoDefinitions_.push_back(std::move(info));
    } else {
        infoDefinitions_.at(found->second) = std::move(info);
    }
    return *this;
}

VcfHeader& VcfHeader::AddSample(std::string sample)
{
    const auto found = sampleLookup_.find(sample);
    if (found == sampleLookup_.cend()) {
        sampleLookup_.insert({sample, samples_.size()});
        samples_.push_back(std::move(sample));
    } else {
        samples_.at(found->second) = std::move(sample);
    }
    return *this;
}

const std::vector<VCF::ContigDefinition>& VcfHeader::ContigDefinitions() const
{
    return contigDefinitions_;
}

const VCF::ContigDefinition& VcfHeader::ContigDefinition(const std::string& id) const
{
    return contigDefinitions_.at(contigLookup_.at(id));
}

VcfHeader& VcfHeader::ContigDefinitions(std::vector<VCF::ContigDefinition> defs)
{
    contigDefinitions_.clear();
    contigLookup_.clear();
    for (auto&& def : defs) {
        AddContigDefinition(std::move(def));
    }
    return *this;
}

const std::string& VcfHeader::FileDate() const
{
    return generalDefinitions_.at(generalLookup_.at("fileDate")).Text();
}

VcfHeader& VcfHeader::FileDate(std::string fileDate)
{
    AddGeneralDefinition({"fileDate", std::move(fileDate)});
    return *this;
}

const std::vector<VCF::FilterDefinition>& VcfHeader::FilterDefinitions() const
{
    return filterDefinitions_;
}

const VCF::FilterDefinition& VcfHeader::FilterDefinition(const std::string& id) const
{
    return filterDefinitions_.at(filterLookup_.at(id));
}

VcfHeader& VcfHeader::FilterDefinitions(std::vector<VCF::FilterDefinition> defs)
{
    filterDefinitions_.clear();
    filterLookup_.clear();
    for (auto&& def : defs) {
        AddFilterDefinition(std::move(def));
    }
    return *this;
}

const std::vector<VCF::FormatDefinition>& VcfHeader::FormatDefinitions() const
{
    return formatDefinitions_;
}

const VCF::FormatDefinition& VcfHeader::FormatDefinition(const std::string& id) const
{
    return formatDefinitions_.at(formatLookup_.at(id));
}

VcfHeader& VcfHeader::FormatDefinitions(std::vector<VCF::FormatDefinition> defs)
{
    formatDefinitions_.clear();
    formatLookup_.clear();
    for (auto&& def : defs) {
        AddFormatDefinition(std::move(def));
    }
    return *this;
}

const std::vector<VCF::GeneralDefinition>& VcfHeader::GeneralDefinitions() const
{
    return generalDefinitions_;
}

const VCF::GeneralDefinition& VcfHeader::GeneralDefinition(const std::string& id) const
{
    return generalDefinitions_.at(generalLookup_.at(id));
}

VcfHeader& VcfHeader::GeneralDefinitions(std::vector<VCF::GeneralDefinition> defs)
{
    generalDefinitions_.clear();
    generalLookup_.clear();
    for (auto&& def : defs) {
        AddGeneralDefinition(std::move(def));
    }
    return *this;
}

const std::vector<VCF::InfoDefinition>& VcfHeader::InfoDefinitions() const
{
    return infoDefinitions_;
}
const VCF::InfoDefinition& VcfHeader::InfoDefinition(const std::string& id) const
{
    return infoDefinitions_.at(infoLookup_.at(id));
}

VcfHeader& VcfHeader::InfoDefinitions(std::vector<VCF::InfoDefinition> defs)
{
    infoDefinitions_.clear();
    infoLookup_.clear();
    for (auto&& def : defs) {
        AddInfoDefinition(std::move(def));
    }
    return *this;
}

size_t VcfHeader::NumLines() const
{
    // +1 for #CHROM line
    return generalDefinitions_.size() + contigDefinitions_.size() + infoDefinitions_.size() +
           filterDefinitions_.size() + formatDefinitions_.size() + 1;
}

const Sample& VcfHeader::SampleAt(size_t index) const { return samples_.at(index); }

size_t VcfHeader::IndexOfSample(const Sample& sample) const { return sampleLookup_.at(sample); }

const std::vector<Sample>& VcfHeader::Samples() const { return samples_; }

VcfHeader& VcfHeader::Samples(std::vector<Sample> names)
{
    samples_.clear();
    sampleLookup_.clear();
    for (auto&& name : names) {
        AddSample(std::move(name));
    }
    return *this;
}

const std::string& VcfHeader::Version() const
{
    return generalDefinitions_.at(generalLookup_.at("fileformat")).Text();
}

VcfHeader& VcfHeader::Version(std::string version)
{
    AddGeneralDefinition({"fileformat", std::move(version)});
    return *this;
}

}  // namespace VCF
}  // namespace PacBio

// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFHEADER_H
#define PBBAM_VCF_VCFHEADER_H

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include <pbbam/vcf/VcfHeaderTypes.h>

namespace PacBio {
namespace VCF {

class VcfHeader
{
public:
    VcfHeader();

    explicit VcfHeader(const std::string& hdrText);

    VcfHeader(const VcfHeader&) = default;
    VcfHeader(VcfHeader&&) = default;
    VcfHeader& operator=(const VcfHeader&) = default;
    VcfHeader& operator=(VcfHeader&&) = default;
    ~VcfHeader() = default;

public:
    // general lines

    size_t NumLines() const;

    const std::string& FileDate() const;
    const std::string& Version() const;

    const std::vector<PacBio::VCF::GeneralDefinition>& GeneralDefinitions() const;
    const PacBio::VCF::GeneralDefinition& GeneralDefinition(const std::string& id) const;

    // ##contig
    const std::vector<PacBio::VCF::ContigDefinition>& ContigDefinitions() const;
    const PacBio::VCF::ContigDefinition& ContigDefinition(const std::string& id) const;

    // INFO

    const std::vector<PacBio::VCF::InfoDefinition>& InfoDefinitions() const;
    const PacBio::VCF::InfoDefinition& InfoDefinition(const std::string& id) const;

    // FILTER

    const std::vector<PacBio::VCF::FilterDefinition>& FilterDefinitions() const;
    const PacBio::VCF::FilterDefinition& FilterDefinition(const std::string& id) const;

    // FORMAT

    const std::vector<PacBio::VCF::FormatDefinition>& FormatDefinitions() const;
    const PacBio::VCF::FormatDefinition& FormatDefinition(const std::string& id) const;

    // samples

    size_t IndexOfSample(const Sample& sample) const;
    const Sample& SampleAt(size_t index) const;
    const std::vector<Sample>& Samples() const;

public:
    // general lines

    VcfHeader& FileDate(std::string fileDate);
    VcfHeader& Version(std::string version);

    VcfHeader& AddGeneralDefinition(PacBio::VCF::GeneralDefinition def);
    VcfHeader& GeneralDefinitions(std::vector<PacBio::VCF::GeneralDefinition> defs);

    // ##contig
    VcfHeader& AddContigDefinition(PacBio::VCF::ContigDefinition def);
    VcfHeader& ContigDefinitions(std::vector<PacBio::VCF::ContigDefinition> defs);

    // INFO

    VcfHeader& AddInfoDefinition(PacBio::VCF::InfoDefinition info);
    VcfHeader& InfoDefinitions(std::vector<PacBio::VCF::InfoDefinition> defs);

    // FILTER

    VcfHeader& AddFilterDefinition(PacBio::VCF::FilterDefinition filter);
    VcfHeader& FilterDefinitions(std::vector<PacBio::VCF::FilterDefinition> defs);

    // FORMAT

    VcfHeader& AddFormatDefinition(PacBio::VCF::FormatDefinition format);
    VcfHeader& FormatDefinitions(std::vector<PacBio::VCF::FormatDefinition> defs);

    // samples

    VcfHeader& AddSample(std::string sample);
    VcfHeader& Samples(std::vector<std::string> names);

private:
    std::vector<PacBio::VCF::GeneralDefinition> generalDefinitions_;
    std::vector<PacBio::VCF::ContigDefinition> contigDefinitions_;
    std::vector<PacBio::VCF::InfoDefinition> infoDefinitions_;
    std::vector<PacBio::VCF::FilterDefinition> filterDefinitions_;
    std::vector<PacBio::VCF::FormatDefinition> formatDefinitions_;
    std::vector<PacBio::VCF::Sample> samples_;

    std::unordered_map<std::string, size_t> generalLookup_;
    std::unordered_map<std::string, size_t> contigLookup_;
    std::unordered_map<std::string, size_t> infoLookup_;
    std::unordered_map<std::string, size_t> filterLookup_;
    std::unordered_map<std::string, size_t> formatLookup_;
    std::unordered_map<std::string, size_t> sampleLookup_;
};

}  // namespace VCF
}  // namespace PacBio

#include "pbbam/vcf/internal/VcfHeader.inl"

#endif  // PBBAM_VCF_VCFHEADER_H

#ifndef PBBAM_VCF_VCFHEADER_H
#define PBBAM_VCF_VCFHEADER_H

#include <pbbam/Config.h>

#include <pbbam/vcf/VcfHeaderTypes.h>

#include <boost/optional.hpp>

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace VCF {

class VcfHeader
{
public:
    VcfHeader();

    explicit VcfHeader(const std::string& hdrText);

public:
    // general lines

    size_t NumLines() const;

    const std::string& FileDate() const;
    const std::string& Version() const;

    const std::vector<VCF::GeneralDefinition>& GeneralDefinitions() const;
    const VCF::GeneralDefinition& GeneralDefinition(const std::string& id) const;

    // ##contig
    const std::vector<VCF::ContigDefinition>& ContigDefinitions() const;
    const VCF::ContigDefinition& ContigDefinition(const std::string& id) const;

    // INFO

    const std::vector<VCF::InfoDefinition>& InfoDefinitions() const;
    const VCF::InfoDefinition& InfoDefinition(const std::string& id) const;

    // FILTER

    const std::vector<VCF::FilterDefinition>& FilterDefinitions() const;
    const VCF::FilterDefinition& FilterDefinition(const std::string& id) const;

    // FORMAT

    const std::vector<VCF::FormatDefinition>& FormatDefinitions() const;
    const VCF::FormatDefinition& FormatDefinition(const std::string& id) const;

    // samples

    size_t IndexOfSample(const Sample& sample) const;
    const Sample& SampleAt(size_t index) const;
    const std::vector<Sample>& Samples() const;

public:
    // general lines

    VcfHeader& FileDate(std::string fileDate);
    VcfHeader& Version(std::string version);

    VcfHeader& AddGeneralDefinition(VCF::GeneralDefinition def);
    VcfHeader& GeneralDefinitions(std::vector<VCF::GeneralDefinition> defs);

    // ##contig
    VcfHeader& AddContigDefinition(VCF::ContigDefinition def);
    VcfHeader& ContigDefinitions(std::vector<VCF::ContigDefinition> defs);

    // INFO

    VcfHeader& AddInfoDefinition(VCF::InfoDefinition info);
    VcfHeader& InfoDefinitions(std::vector<VCF::InfoDefinition> defs);

    // FILTER

    VcfHeader& AddFilterDefinition(VCF::FilterDefinition filter);
    VcfHeader& FilterDefinitions(std::vector<VCF::FilterDefinition> defs);

    // FORMAT

    VcfHeader& AddFormatDefinition(VCF::FormatDefinition format);
    VcfHeader& FormatDefinitions(std::vector<VCF::FormatDefinition> defs);

    // samples

    VcfHeader& AddSample(std::string sample);
    VcfHeader& Samples(std::vector<std::string> names);

private:
    std::vector<VCF::GeneralDefinition> generalDefinitions_;
    std::vector<VCF::ContigDefinition> contigDefinitions_;
    std::vector<VCF::InfoDefinition> infoDefinitions_;
    std::vector<VCF::FilterDefinition> filterDefinitions_;
    std::vector<VCF::FormatDefinition> formatDefinitions_;
    std::vector<VCF::Sample> samples_;

    std::unordered_map<std::string, size_t> generalLookup_;
    std::unordered_map<std::string, size_t> contigLookup_;
    std::unordered_map<std::string, size_t> infoLookup_;
    std::unordered_map<std::string, size_t> filterLookup_;
    std::unordered_map<std::string, size_t> formatLookup_;
    std::unordered_map<std::string, size_t> sampleLookup_;
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFHEADER_H

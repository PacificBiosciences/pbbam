// Author: Derek Barnett

#ifndef PBBAM_VCF_VARIANT_H
#define PBBAM_VCF_VARIANT_H

#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include <pbbam/Position.h>

namespace PacBio {
namespace VCF {

struct InfoField
{
    std::string id;
    boost::optional<std::string> value;
    boost::optional<std::vector<std::string>> values;
};

struct GenotypeData
{
    boost::optional<std::string> value;
    boost::optional<std::vector<std::string>> values;
};

struct GenotypeField
{
    std::vector<GenotypeData> data;
};

class VcfVariant
{
public:
    VcfVariant();

    explicit VcfVariant(const std::string& text);

    VcfVariant(std::string id, std::string chrom, PacBio::BAM::Position pos, std::string refAllele,
               std::string altAllele);

    VcfVariant(const VcfVariant&) = default;
    VcfVariant(VcfVariant&&) = default;
    VcfVariant& operator=(const VcfVariant&) = default;
    VcfVariant& operator=(VcfVariant&&) = default;
    ~VcfVariant() = default;

public:
    // core fields

    const std::string& Chrom() const;
    VcfVariant& Chrom(std::string chrom);

    PacBio::BAM::Position Position() const;
    VcfVariant& Position(PacBio::BAM::Position pos);

    const std::string& Id() const;
    VcfVariant& Id(std::string id);

    const std::string& RefAllele() const;
    VcfVariant& RefAllele(std::string refAllele);

    const std::string& AltAllele() const;
    VcfVariant& AltAllele(std::string altAllele);

    float Quality() const;
    VcfVariant& Quality(float qual);

    const std::string& Filter() const;
    VcfVariant& Filter(std::string filter);

    // convenience methods
    bool IsDeletion() const;
    bool IsInsertion() const;
    bool IsQualityMissing() const;
    bool IsSnp() const;

public:
    // info fields

    VcfVariant& AddInfoField(InfoField field);
    VcfVariant& RemoveInfoField(const std::string& id);

    const std::vector<InfoField>& InfoFields() const;
    VcfVariant& InfoFields(std::vector<InfoField> fields);

    bool HasInfoField(const std::string& id) const;

    const boost::optional<std::string> InfoValue(const std::string& id) const;
    VcfVariant& InfoValue(const std::string& id, boost::optional<std::string> value);

    const boost::optional<std::vector<std::string>> InfoValues(const std::string& id) const;
    VcfVariant& InfoValues(const std::string& id, boost::optional<std::vector<std::string>> values);

public:
    // sample genotypes

    // NOTE: if you want to look up by sample name, get the index from header

    std::vector<std::string> GenotypeIds() const;
    VcfVariant& GenotypeIds(std::vector<std::string> ids);

    std::vector<GenotypeField> Genotypes() const;
    VcfVariant& Genotypes(std::vector<GenotypeField> genotypes);

    const boost::optional<std::string>& GenotypeValue(const size_t sampleIndex,
                                                      const std::string& id) const;
    VcfVariant& GenotypeValue(const size_t sampleIndex, const std::string& id,
                              boost::optional<std::string> value);

    const boost::optional<std::vector<std::string>>& GenotypeValues(const size_t sampleIndex,
                                                                    const std::string& id) const;
    VcfVariant& GenotypeValues(const size_t sampleIndex, const std::string& id,
                               boost::optional<std::vector<std::string>> values);

    bool IsSampleHeterozygous(const size_t sampleIndex) const;
    bool IsSamplePhased(const size_t sampleIndex) const;

private:
    // FIXED data
    std::string chrom_;
    PacBio::BAM::Position pos_;
    std::string id_;
    std::string refAllele_;
    std::string altAllele_;  // multiple? KISS, only add if needed
    float qual_;
    std::string filter_;

    // INFO data
    std::vector<InfoField> infoFields_;
    std::unordered_map<std::string, size_t> infoLookup_;

    // SAMPLE GENOTYPE data
    std::vector<std::string> format_;  // order matches FORMAT string
    std::unordered_map<std::string, size_t>
        genotypeDataLookup_;                      // genotype ID -> genotypeField.data index
    std::vector<GenotypeField> sampleGenotypes_;  // index matches sample order
};

}  // namespace VCF
}  // namespace PacBio

#include "pbbam/vcf/internal/VcfVariant.inl"

#endif  // PBBAM_VCF_VARIANT_H

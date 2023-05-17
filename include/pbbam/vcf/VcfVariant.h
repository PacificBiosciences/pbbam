#ifndef PBBAM_VCF_VARIANT_H
#define PBBAM_VCF_VARIANT_H

#include <pbbam/Config.h>

#include <pbbam/vcf/VcfHeaderTypes.h>

#include <pbcopper/data/Position.h>

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace PacBio {
namespace VCF {

struct InfoField
{
    std::string id;
    std::optional<std::string> value;
    std::optional<std::vector<std::string>> values;
};

struct GenotypeData
{
    std::optional<std::string> value;
    std::optional<std::vector<std::string>> values;
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

    VcfVariant(std::string id, std::string chrom, Data::Position pos, std::string refAllele,
               std::string altAllele);

public:
    // core fields

    const std::string& Chrom() const;
    VcfVariant& Chrom(std::string chrom);

    Data::Position Position() const;
    VcfVariant& Position(Data::Position pos);

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

    const std::optional<std::string> InfoValue(const std::string& id) const;
    VcfVariant& InfoValue(const std::string& id, std::optional<std::string> value);

    const std::optional<std::vector<std::string>> InfoValues(const std::string& id) const;
    VcfVariant& InfoValues(const std::string& id, std::optional<std::vector<std::string>> values);

public:
    // sample genotypes

    // NOTE: if you want to look up by sample name, get the index from header

    std::vector<std::string> GenotypeIds() const;
    VcfVariant& GenotypeIds(std::vector<std::string> ids);

    std::vector<GenotypeField> Genotypes() const;
    VcfVariant& Genotypes(std::vector<GenotypeField> genotypes);

    const std::optional<std::string>& GenotypeValue(std::size_t sampleIndex,
                                                    const std::string& id) const;
    VcfVariant& GenotypeValue(std::size_t sampleIndex, const std::string& id,
                              std::optional<std::string> value);

    const std::optional<std::vector<std::string>>& GenotypeValues(std::size_t sampleIndex,
                                                                  const std::string& id) const;
    VcfVariant& GenotypeValues(std::size_t sampleIndex, const std::string& id,
                               std::optional<std::vector<std::string>> values);

    bool IsSampleHeterozygous(std::size_t sampleIndex) const;
    bool IsSamplePhased(std::size_t sampleIndex) const;

private:
    // FIXED data
    std::string chrom_;
    Data::Position pos_;
    std::string id_;
    std::string refAllele_;
    std::string altAllele_;  // multiple? KISS, only add if needed
    float qual_;
    std::string filter_;

    // INFO data
    std::vector<InfoField> infoFields_;
    std::unordered_map<std::string, std::size_t> infoLookup_;

    // SAMPLE GENOTYPE data
    std::vector<std::string> format_;  // order matches FORMAT string
    std::unordered_map<std::string, std::size_t>
        genotypeDataLookup_;                      // genotype ID -> genotypeField.data index
    std::vector<GenotypeField> sampleGenotypes_;  // index matches sample order
};

}  // namespace VCF
}  // namespace PacBio

// #include "pbbam/vcf/internal/VcfVariant.inl"

#endif  // PBBAM_VCF_VARIANT_H

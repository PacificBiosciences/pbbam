#ifndef PBBAM_VCF_VCFVcfVariant_INL
#define PBBAM_VCF_VCFVcfVariant_INL

#include <pbbam/vcf/VcfVariant.h>

#include <cmath>

#include <pbbam/StringUtilities.h>

namespace PacBio {
namespace VCF {

inline  VcfVariant::VcfVariant()
    : pos_{PacBio::BAM::UnmappedPosition}
    , qual_{NAN}
    , filter_{"PASS"}
{
}

inline VcfVariant::VcfVariant(
        std::string id, std::string chrom, PacBio::BAM::Position pos,
        std::string refAllele, std::string altAllele)
    : chrom_{std::move(chrom)}
    , pos_{pos}
    , id_{std::move(id)}
    , refAllele_{std::move(refAllele)}
    , altAllele_{std::move(altAllele)}
    , qual_{NAN}
    , filter_{"PASS"}
{
}

inline VcfVariant& VcfVariant::AddInfoField(InfoField field)
{
    const auto found = infoLookup_.find(field.id);
    if (found == infoLookup_.cend()) {
        infoLookup_.insert({field.id, infoFields_.size()});
        infoFields_.push_back(std::move(field));
    }
    else
        infoFields_.at(found->second) = std::move(field);
    return *this;
}

inline const std::string& VcfVariant::AltAllele() const { return altAllele_; }

inline VcfVariant& VcfVariant::AltAllele(std::string altAllele)
{
    altAllele_ = std::move(altAllele);
    return *this;
}

inline const std::string& VcfVariant::Chrom() const { return chrom_; }

inline VcfVariant& VcfVariant::Chrom(std::string chrom)
{
    chrom_ = std::move(chrom);
    return *this;
}

inline const std::string& VcfVariant::Filter() const { return filter_; }

inline VcfVariant& VcfVariant::Filter(std::string filter)
{
    filter_ = std::move(filter);
    return *this;
}

inline std::vector<std::string> VcfVariant::GenotypeIds() const
{
    return format_;
}

inline VcfVariant& VcfVariant::GenotypeIds(std::vector<std::string> ids)
{
    genotypeDataLookup_.clear();

    format_ = std::move(ids);
    for (size_t i = 0; i < format_.size(); ++i)
        genotypeDataLookup_.insert({format_.at(i), i});
    return *this;
}

inline std::vector<GenotypeField> VcfVariant::Genotypes() const
{
    return sampleGenotypes_;
}

inline VcfVariant& VcfVariant::Genotypes(std::vector<GenotypeField> genotypes)
{
    sampleGenotypes_ = std::move(genotypes);
    return *this;
}

inline const boost::optional<std::string>& VcfVariant::GenotypeValue(
        const size_t sampleIndex, const std::string& id) const
{
    const auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    const auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    return genotypeData.value;
}

inline VcfVariant& VcfVariant::GenotypeValue(
        const size_t sampleIndex, const std::string& id, boost::optional<std::string> value)
{
    auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    genotypeData.value = std::move(value);
    return *this;
}

inline const boost::optional<std::vector<std::string>>& VcfVariant::GenotypeValues(
        const size_t sampleIndex, const std::string& id) const
{
    const auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    const auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    return genotypeData.values;
}

inline VcfVariant& VcfVariant::GenotypeValues(
        const size_t sampleIndex, const std::string& id,
        boost::optional<std::vector<std::string>> values)
{
    auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    genotypeData.values = std::move(values);
    return *this;
}

inline bool VcfVariant::HasInfoField(const std::string& id) const
{
    const auto found = infoLookup_.find(id);
    return found != infoLookup_.cend();
}

inline const std::string& VcfVariant::Id() const { return id_; }

inline VcfVariant& VcfVariant::Id(std::string id)
{
    id_ = std::move(id);
    return *this;
}

inline const std::vector<InfoField>& VcfVariant::InfoFields() const
{
    return infoFields_;
}

inline VcfVariant& VcfVariant::InfoFields(std::vector<InfoField> fields)
{
    infoFields_.clear();
    infoLookup_.clear();
    for (auto&& field : fields)
        AddInfoField(std::move(field));
    return *this;
}

inline const boost::optional<std::string> VcfVariant::InfoValue(const std::string& id) const
{
    return infoFields_.at(infoLookup_.at(id)).value;
}

inline VcfVariant& VcfVariant::InfoValue(const std::string& id, boost::optional<std::string> value)
{
    infoFields_.at(infoLookup_.at(id)).value = std::move(value);
    return *this;
}

inline const boost::optional<std::vector<std::string>> VcfVariant::InfoValues(const std::string& id) const
{
    return infoFields_.at(infoLookup_.at(id)).values;
}

inline VcfVariant& VcfVariant::InfoValues(const std::string& id, boost::optional<std::vector<std::string>> values)
{
    infoFields_.at(infoLookup_.at(id)).values = std::move(values);
    return *this;
}

inline bool VcfVariant::IsDeletion() const
{
    return refAllele_.size() > altAllele_.size();
}

inline bool VcfVariant::IsInsertion() const
{
    return refAllele_.size() < altAllele_.size();
}

inline bool VcfVariant::IsQualityMissing() const
{
    return std::isnan(qual_);
}

inline bool VcfVariant::IsSampleHeterozygous(const size_t sampleIndex) const
{
    const auto data = GenotypeValue(sampleIndex, "GT");
    auto fields = PacBio::BAM::Split(data.get(), '/');
    if (fields.size() == 1)
        fields = PacBio::BAM::Split(data.get(), '|');

    if (fields.size() != 2)
        throw std::runtime_error{"VCF format error: malformatted GT field: " + data.get()};

    return fields.at(0) != fields.at(1);
}

inline bool VcfVariant::IsSamplePhased(const size_t sampleIndex) const
{
    const auto data = GenotypeValue(sampleIndex, "GT");
    const auto phaseFound = data.get().find('|') != std::string::npos;
    if (phaseFound) assert(data.get().find('/') == std::string::npos);
    return phaseFound;
}

inline bool VcfVariant::IsSnp() const
{
    return refAllele_.size() == 1 &&
           altAllele_.size() == 1 &&
           refAllele_[0] != altAllele_[0];
}

inline PacBio::BAM::Position VcfVariant::Position() const { return pos_; }

inline VcfVariant& VcfVariant::Position(PacBio::BAM::Position pos)
{
    pos_ = pos;
    return *this;
}

inline float VcfVariant::Quality() const { return qual_; }

inline VcfVariant& VcfVariant::Quality(float qual)
{
    qual_ = qual;
    return *this;
}

inline const std::string& VcfVariant::RefAllele() const { return refAllele_; }

inline VcfVariant& VcfVariant::RefAllele(std::string refAllele)
{
    refAllele_ = std::move(refAllele);
    return *this;
}

inline VcfVariant& VcfVariant::RemoveInfoField(const std::string& id)
{
    const auto found = infoLookup_.find(id);
    if (found == infoLookup_.cend())
        return *this;

    const auto currentFields = InfoFields();

    infoFields_.clear();
    infoLookup_.clear();

    for (auto&& field : currentFields)
    {
        if (field.id != id)
            AddInfoField(std::move(field));
    }

    return *this;
}

} // namespace VCF
} // namespace PacBio

#endif // PBBAM_VCF_VCFVcfVariant_INL

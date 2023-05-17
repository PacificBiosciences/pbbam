#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfVariant.h>

#include <pbbam/StringUtilities.h>
#include <pbbam/vcf/VcfFormat.h>
#include "VcfFormatException.h"

#include <type_traits>

#include <cassert>
#include <cmath>

namespace PacBio {
namespace VCF {

VcfVariant::VcfVariant(const std::string& text) { *this = VcfFormat::ParsedVariant(text); }

VcfVariant::VcfVariant() : pos_{Data::UNMAPPED_POSITION}, qual_{NAN}, filter_{"PASS"} {}

VcfVariant::VcfVariant(std::string id, std::string chrom, Data::Position pos, std::string refAllele,
                       std::string altAllele)
    : chrom_{std::move(chrom)}
    , pos_{pos}
    , id_{std::move(id)}
    , refAllele_{std::move(refAllele)}
    , altAllele_{std::move(altAllele)}
    , qual_{NAN}
    , filter_{"PASS"}
{}

VcfVariant& VcfVariant::AddInfoField(InfoField field)
{
    const auto found = infoLookup_.find(field.id);
    if (found == infoLookup_.cend()) {
        infoLookup_.insert({field.id, infoFields_.size()});
        infoFields_.push_back(std::move(field));
    } else {
        infoFields_.at(found->second) = std::move(field);
    }
    return *this;
}

const std::string& VcfVariant::AltAllele() const { return altAllele_; }

VcfVariant& VcfVariant::AltAllele(std::string altAllele)
{
    altAllele_ = std::move(altAllele);
    return *this;
}

const std::string& VcfVariant::Chrom() const { return chrom_; }

VcfVariant& VcfVariant::Chrom(std::string chrom)
{
    chrom_ = std::move(chrom);
    return *this;
}

const std::string& VcfVariant::Filter() const { return filter_; }

VcfVariant& VcfVariant::Filter(std::string filter)
{
    filter_ = std::move(filter);
    return *this;
}

std::vector<std::string> VcfVariant::GenotypeIds() const { return format_; }

VcfVariant& VcfVariant::GenotypeIds(std::vector<std::string> ids)
{
    genotypeDataLookup_.clear();

    format_ = std::move(ids);
    for (std::size_t i = 0; i < format_.size(); ++i) {
        genotypeDataLookup_.emplace(format_.at(i), i);
    }
    return *this;
}

std::vector<GenotypeField> VcfVariant::Genotypes() const { return sampleGenotypes_; }

VcfVariant& VcfVariant::Genotypes(std::vector<GenotypeField> genotypes)
{
    sampleGenotypes_ = std::move(genotypes);
    return *this;
}

const std::optional<std::string>& VcfVariant::GenotypeValue(const std::size_t sampleIndex,
                                                            const std::string& id) const
{
    const auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    const auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    return genotypeData.value;
}

VcfVariant& VcfVariant::GenotypeValue(const std::size_t sampleIndex, const std::string& id,
                                      std::optional<std::string> value)
{
    auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    genotypeData.value = std::move(value);
    return *this;
}

const std::optional<std::vector<std::string>>& VcfVariant::GenotypeValues(
    const std::size_t sampleIndex, const std::string& id) const
{
    const auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    const auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    return genotypeData.values;
}

VcfVariant& VcfVariant::GenotypeValues(const std::size_t sampleIndex, const std::string& id,
                                       std::optional<std::vector<std::string>> values)
{
    auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    genotypeData.values = std::move(values);
    return *this;
}

bool VcfVariant::HasInfoField(const std::string& id) const
{
    const auto found = infoLookup_.find(id);
    return found != infoLookup_.cend();
}

const std::string& VcfVariant::Id() const { return id_; }

VcfVariant& VcfVariant::Id(std::string id)
{
    id_ = std::move(id);
    return *this;
}

const std::vector<InfoField>& VcfVariant::InfoFields() const { return infoFields_; }

VcfVariant& VcfVariant::InfoFields(std::vector<InfoField> fields)
{
    infoFields_.clear();
    infoLookup_.clear();
    for (auto&& field : fields) {
        AddInfoField(std::move(field));
    }
    return *this;
}

const std::optional<std::string> VcfVariant::InfoValue(const std::string& id) const
{
    return infoFields_.at(infoLookup_.at(id)).value;
}

VcfVariant& VcfVariant::InfoValue(const std::string& id, std::optional<std::string> value)
{
    infoFields_.at(infoLookup_.at(id)).value = std::move(value);
    return *this;
}

const std::optional<std::vector<std::string>> VcfVariant::InfoValues(const std::string& id) const
{
    return infoFields_.at(infoLookup_.at(id)).values;
}

VcfVariant& VcfVariant::InfoValues(const std::string& id,
                                   std::optional<std::vector<std::string>> values)
{
    infoFields_.at(infoLookup_.at(id)).values = std::move(values);
    return *this;
}

bool VcfVariant::IsDeletion() const { return refAllele_.size() > altAllele_.size(); }

bool VcfVariant::IsInsertion() const { return refAllele_.size() < altAllele_.size(); }

bool VcfVariant::IsQualityMissing() const { return std::isnan(qual_); }

bool VcfVariant::IsSampleHeterozygous(const std::size_t sampleIndex) const
{
    const auto data = GenotypeValue(sampleIndex, "GT");
    auto fields = BAM::Split(*data, '/');
    if (fields.size() == 1) {
        fields = BAM::Split(*data, '|');
    }

    if (fields.size() == 2) {
        return fields.at(0) != fields.at(1);
    } else {
        throw VcfFormatException{"malformed GT field: " + *data};
    }
}

bool VcfVariant::IsSamplePhased(const std::size_t sampleIndex) const
{
    const auto data = GenotypeValue(sampleIndex, "GT");
    const auto phaseFound = data->find('|') != std::string::npos;
    if (phaseFound) {
        assert(data->find('/') == std::string::npos);
    }
    return phaseFound;
}

bool VcfVariant::IsSnp() const
{
    return refAllele_.size() == 1 && altAllele_.size() == 1 && refAllele_[0] != altAllele_[0];
}

Data::Position VcfVariant::Position() const { return pos_; }

VcfVariant& VcfVariant::Position(Data::Position pos)
{
    pos_ = pos;
    return *this;
}

float VcfVariant::Quality() const { return qual_; }

VcfVariant& VcfVariant::Quality(float qual)
{
    qual_ = qual;
    return *this;
}

const std::string& VcfVariant::RefAllele() const { return refAllele_; }

VcfVariant& VcfVariant::RefAllele(std::string refAllele)
{
    refAllele_ = std::move(refAllele);
    return *this;
}

VcfVariant& VcfVariant::RemoveInfoField(const std::string& id)
{
    const auto found = infoLookup_.find(id);
    if (found == infoLookup_.cend()) {
        return *this;
    }

    const auto currentFields = InfoFields();

    infoFields_.clear();
    infoLookup_.clear();

    for (auto&& field : currentFields) {
        if (field.id != id) {
            AddInfoField(std::move(field));
        }
    }

    return *this;
}

}  // namespace VCF
}  // namespace PacBio

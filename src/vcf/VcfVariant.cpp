// Author: Derek Barnett

#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfVariant.h>

#include <cassert>
#include <cmath>
#include <type_traits>

#include <pbbam/StringUtilities.h>
#include <pbbam/vcf/VcfFormat.h>

namespace PacBio {
namespace VCF {

static_assert(std::is_copy_constructible<VcfVariant>::value,
              "VcfVariant(const VcfVariant&) is not = default");
static_assert(std::is_copy_assignable<VcfVariant>::value,
              "VcfVariant& operator=(const VcfVariant&) is not = default");

static_assert(std::is_nothrow_move_constructible<VcfVariant>::value,
              "VcfVariant(VcfVariant&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<VcfVariant>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

VcfVariant::VcfVariant(const std::string& text) { *this = VcfFormat::ParsedVariant(text); }

VcfVariant::VcfVariant() : pos_{PacBio::BAM::UnmappedPosition}, qual_{NAN}, filter_{"PASS"} {}

VcfVariant::VcfVariant(std::string id, std::string chrom, PacBio::BAM::Position pos,
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

VcfVariant& VcfVariant::AddInfoField(InfoField field)
{
    const auto found = infoLookup_.find(field.id);
    if (found == infoLookup_.cend()) {
        infoLookup_.insert({field.id, infoFields_.size()});
        infoFields_.push_back(std::move(field));
    } else
        infoFields_.at(found->second) = std::move(field);
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
    for (size_t i = 0; i < format_.size(); ++i)
        genotypeDataLookup_.insert({format_.at(i), i});
    return *this;
}

std::vector<GenotypeField> VcfVariant::Genotypes() const { return sampleGenotypes_; }

VcfVariant& VcfVariant::Genotypes(std::vector<GenotypeField> genotypes)
{
    sampleGenotypes_ = std::move(genotypes);
    return *this;
}

const boost::optional<std::string>& VcfVariant::GenotypeValue(const size_t sampleIndex,
                                                              const std::string& id) const
{
    const auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    const auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    return genotypeData.value;
}

VcfVariant& VcfVariant::GenotypeValue(const size_t sampleIndex, const std::string& id,
                                      boost::optional<std::string> value)
{
    auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    genotypeData.value = std::move(value);
    return *this;
}

const boost::optional<std::vector<std::string>>& VcfVariant::GenotypeValues(
    const size_t sampleIndex, const std::string& id) const
{
    const auto& genotypeField = sampleGenotypes_.at(sampleIndex);
    const auto genotypeDataIndex = genotypeDataLookup_.at(id);
    const auto& genotypeData = genotypeField.data.at(genotypeDataIndex);
    return genotypeData.values;
}

VcfVariant& VcfVariant::GenotypeValues(const size_t sampleIndex, const std::string& id,
                                       boost::optional<std::vector<std::string>> values)
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
    for (auto&& field : fields)
        AddInfoField(std::move(field));
    return *this;
}

const boost::optional<std::string> VcfVariant::InfoValue(const std::string& id) const
{
    return infoFields_.at(infoLookup_.at(id)).value;
}

VcfVariant& VcfVariant::InfoValue(const std::string& id, boost::optional<std::string> value)
{
    infoFields_.at(infoLookup_.at(id)).value = std::move(value);
    return *this;
}

const boost::optional<std::vector<std::string>> VcfVariant::InfoValues(const std::string& id) const
{
    return infoFields_.at(infoLookup_.at(id)).values;
}

VcfVariant& VcfVariant::InfoValues(const std::string& id,
                                   boost::optional<std::vector<std::string>> values)
{
    infoFields_.at(infoLookup_.at(id)).values = std::move(values);
    return *this;
}

bool VcfVariant::IsDeletion() const { return refAllele_.size() > altAllele_.size(); }

bool VcfVariant::IsInsertion() const { return refAllele_.size() < altAllele_.size(); }

bool VcfVariant::IsQualityMissing() const { return std::isnan(qual_); }

bool VcfVariant::IsSampleHeterozygous(const size_t sampleIndex) const
{
    const auto data = GenotypeValue(sampleIndex, "GT");
    auto fields = PacBio::BAM::Split(data.get(), '/');
    if (fields.size() == 1) fields = PacBio::BAM::Split(data.get(), '|');

    if (fields.size() != 2)
        throw std::runtime_error{"VcfFormat: malformatted GT field: " + data.get()};

    return fields.at(0) != fields.at(1);
}

bool VcfVariant::IsSamplePhased(const size_t sampleIndex) const
{
    const auto data = GenotypeValue(sampleIndex, "GT");
    const auto phaseFound = data.get().find('|') != std::string::npos;
    if (phaseFound) assert(data.get().find('/') == std::string::npos);
    return phaseFound;
}

bool VcfVariant::IsSnp() const
{
    return refAllele_.size() == 1 && altAllele_.size() == 1 && refAllele_[0] != altAllele_[0];
}

PacBio::BAM::Position VcfVariant::Position() const { return pos_; }

VcfVariant& VcfVariant::Position(PacBio::BAM::Position pos)
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
    if (found == infoLookup_.cend()) return *this;

    const auto currentFields = InfoFields();

    infoFields_.clear();
    infoLookup_.clear();

    for (auto&& field : currentFields) {
        if (field.id != id) AddInfoField(std::move(field));
    }

    return *this;
}

}  // namespace VCF
}  // namespace PacBio

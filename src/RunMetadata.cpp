// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/RunMetadata.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <tuple>

#include <boost/algorithm/string/predicate.hpp>

#include "RunMetadataParser.h"

namespace PacBio {
namespace BAM {

// -----------------------
// CollectionElementBase
// -----------------------

CollectionElementBase::CollectionElementBase() = default;

CollectionElementBase::CollectionElementBase(std::map<std::string, std::string> data)
    : data_{std::move(data)}
{
}

const std::string& CollectionElementBase::Get(const std::string& name) const
{
    return data_.at(name);
}

std::string& CollectionElementBase::Get(const std::string& name) { return data_.at(name); }

bool CollectionElementBase::Has(const std::string& name) const
{
    return (data_.find(name) != data_.cend());
}

CollectionElementBase& CollectionElementBase::Set(const std::string& name, const std::string& value)
{
    data_[name] = value;
    return *this;
}

// ----------------------
// AutomationParameters
// ----------------------

AutomationParameters::AutomationParameters() = default;

AutomationParameters::AutomationParameters(std::map<std::string, std::string> data)
    : CollectionElementBase{std::move(data)}
{
}

int32_t AutomationParameters::CellNFCIndex() const { return std::stoi(Get(Element::CellNFCIndex)); }

AutomationParameters& AutomationParameters::CellNFCIndex(int32_t i)
{
    Set(Element::CellNFCIndex, std::to_string(i));
    return *this;
}

bool AutomationParameters::HasCellNFCIndex() const { return Has(Element::CellNFCIndex); }

int32_t AutomationParameters::CollectionNumber() const
{
    return std::stoi(Get(Element::CollectionNumber));
}

AutomationParameters& AutomationParameters::CollectionNumber(int32_t i)
{
    Set(Element::CollectionNumber, std::to_string(i));
    return *this;
}

bool AutomationParameters::HasCollectionNumber() const { return Has(Element::CellNFCIndex); }

double AutomationParameters::Exposure() const { return std::stod(Get(Element::Exposure)); }

AutomationParameters& AutomationParameters::Exposure(double d)
{
    Set(Element::Exposure, std::to_string(d));
    return *this;
}

bool AutomationParameters::HasExposure() const { return Has(Element::Exposure); }

bool AutomationParameters::ExtendFirst() const
{
    return boost::iequals(Get(Element::ExtendFirst), "True");
}

AutomationParameters& AutomationParameters::ExtendFirst(bool ok)
{
    Set(Element::ExtendFirst, (ok ? "True" : "False"));
    return *this;
}

bool AutomationParameters::HasExtendFirst() const { return Has(Element::ExtendFirst); }

double AutomationParameters::ExtensionTime() const
{
    return std::stod(Get(Element::ExtensionTime));
}

AutomationParameters& AutomationParameters::ExtensionTime(double d)
{
    Set(Element::ExtensionTime, std::to_string(d));
    return *this;
}

bool AutomationParameters::HasExtensionTime() const { return Has(Element::ExtensionTime); }

int32_t AutomationParameters::ExtraIMWashes() const
{
    return std::stoi(Get(Element::ExtraIMWashes));
}

AutomationParameters& AutomationParameters::ExtraIMWashes(int32_t i)
{
    Set(Element::ExtraIMWashes, std::to_string(i));
    return *this;
}

bool AutomationParameters::HasExtraIMWashes() const { return Has(Element::ExtraIMWashes); }

bool AutomationParameters::HasN2Switch() const
{
    return boost::iequals(Get(Element::HasN2Switch), "True");
}

AutomationParameters& AutomationParameters::HasN2Switch(bool ok)
{
    Set(Element::HasN2Switch, (ok ? "True" : "False"));
    return *this;
}

bool AutomationParameters::HasHasN2Switch() const { return Has(Element::HasN2Switch); }

std::string AutomationParameters::HQRFMethod() const { return Get(Element::HQRFMethod); }

AutomationParameters& AutomationParameters::HQRFMethod(std::string s)
{
    Set(Element::HQRFMethod, s);
    return *this;
}

bool AutomationParameters::HasHQRFMethod() const { return Has(Element::HQRFMethod); }

double AutomationParameters::ImmobilizationTime() const
{
    return std::stod(Get(Element::ImmobilizationTime));
}

AutomationParameters& AutomationParameters::ImmobilizationTime(double d)
{
    Set(Element::ImmobilizationTime, std::to_string(d));
    return *this;
}

bool AutomationParameters::HasImmobilizationTime() const
{
    return Has(Element::ImmobilizationTime);
}

int32_t AutomationParameters::InsertSize() const { return std::stoi(Get(Element::InsertSize)); }

AutomationParameters& AutomationParameters::InsertSize(int32_t i)
{
    Set(Element::InsertSize, std::to_string(i));
    return *this;
}
bool AutomationParameters::HasInsertSize() const { return Has(Element::InsertSize); }

double AutomationParameters::MovieLength() const { return std::stod(Get(Element::MovieLength)); }

AutomationParameters& AutomationParameters::MovieLength(double d)
{
    Set(Element::ImmobilizationTime, std::to_string(d));
    return *this;
}

bool AutomationParameters::HasMovieLength() const { return Has(Element::MovieLength); }

bool AutomationParameters::PCDinPlate() const
{
    return boost::iequals(Get(Element::PCDinPlate), "True");
}

AutomationParameters& AutomationParameters::PCDinPlate(bool ok)
{
    Set(Element::PCDinPlate, (ok ? "True" : "False"));
    return *this;
}

bool AutomationParameters::HasPCDinPlate() const { return Has(Element::PCDinPlate); }

bool AutomationParameters::PreExtensionWorkflow() const
{
    return boost::iequals(Get(Element::PreExtensionWorkflow), "True");
}

AutomationParameters& AutomationParameters::PreExtensionWorkflow(bool ok)
{
    Set(Element::PreExtensionWorkflow, (ok ? "True" : "False"));
    return *this;
}

bool AutomationParameters::HasPreExtensionWorkflow() const
{
    return Has(Element::PreExtensionWorkflow);
}

double AutomationParameters::SNRCut() const { return std::stod(Get(Element::SNRCut)); }

AutomationParameters& AutomationParameters::SNRCut(double d)
{
    Set(Element::SNRCut, std::to_string(d));
    return *this;
}

bool AutomationParameters::HasSNRCut() const { return Has(Element::SNRCut); }

int32_t AutomationParameters::TipSearchMaxDuration() const
{
    return std::stoi(Get(Element::TipSearchMaxDuration));
}

AutomationParameters& AutomationParameters::TipSearchMaxDuration(int32_t i)
{
    Set(Element::TipSearchMaxDuration, std::to_string(i));
    return *this;
}

bool AutomationParameters::HasTipSearchMaxDuration() const
{
    return Has(Element::TipSearchMaxDuration);
}

bool AutomationParameters::UseStageHotStart() const
{
    return boost::iequals(Get(Element::UseStageHotStart), "True");
}

AutomationParameters& AutomationParameters::UseStageHotStart(bool ok)
{
    Set(Element::UseStageHotStart, (ok ? "True" : "False"));
    return *this;
}

bool AutomationParameters::HasUseStageHotStart() const { return Has(Element::UseStageHotStart); }

std::string AutomationParameters::GetParameter(const std::string& param) const
{
    return Get(param);
}

AutomationParameters& AutomationParameters::SetParameter(const std::string& name,
                                                         const std::string& value)
{
    Set(name, value);
    return *this;
}

bool AutomationParameters::HasParameter(const std::string& param) const { return Has(param); }

// ----------------------
// BindingKit
// ----------------------

BindingKit::BindingKit() = default;

BindingKit::BindingKit(std::map<std::string, std::string> data)
    : CollectionElementBase{std::move(data)}
{
}

std::string BindingKit::PartNumber() const { return Get(Element::PartNumber); }

BindingKit& BindingKit::PartNumber(std::string s)
{
    Set(Element::PartNumber, s);
    return *this;
}

bool BindingKit::HasPartNumber() { return Has(Element::PartNumber); }

// ----------------------
// ControlKit
// ----------------------

ControlKit::ControlKit() = default;

ControlKit::ControlKit(std::map<std::string, std::string> data)
    : CollectionElementBase{std::move(data)}
{
}

std::string ControlKit::PartNumber() const { return Get(Element::PartNumber); }

ControlKit& ControlKit::PartNumber(std::string s)
{
    Set(Element::PartNumber, s);
    return *this;
}

bool ControlKit::HasPartNumber() const { return Has(Element::PartNumber); }

std::string ControlKit::LeftAdapter() const { return Get(Element::LeftAdapter); }

ControlKit& ControlKit::LeftAdapter(std::string s)
{
    Set(Element::LeftAdapter, s);
    return *this;
}

bool ControlKit::HasLeftAdapter() const { return Has(Element::LeftAdapter); }

std::string ControlKit::RightAdapter() const { return Get(Element::RightAdapter); }

ControlKit& ControlKit::RightAdapter(std::string s)
{
    Set(Element::RightAdapter, s);
    return *this;
}

bool ControlKit::HasRightAdapter() const { return Has(Element::RightAdapter); }

std::string ControlKit::Sequence() const { return Get(Element::Sequence); }

ControlKit& ControlKit::Sequence(std::string s)
{
    Set(Element::Sequence, s);
    return *this;
}

bool ControlKit::HasSequence() const { return Has(Element::Sequence); }

// ----------------------
// SequencingKitPlate
// ----------------------

SequencingKitPlate::SequencingKitPlate() = default;

SequencingKitPlate::SequencingKitPlate(std::map<std::string, std::string> data)
    : CollectionElementBase{std::move(data)}
{
}

std::string SequencingKitPlate::PartNumber() const { return Get(Element::PartNumber); }

SequencingKitPlate& SequencingKitPlate::PartNumber(std::string s)
{
    Set(Element::PartNumber, s);
    return *this;
}

bool SequencingKitPlate::HasPartNumber() const { return Has(Element::PartNumber); }

// ----------------------
// TemplatePrepKit
// ----------------------

TemplatePrepKit::TemplatePrepKit() = default;

TemplatePrepKit::TemplatePrepKit(std::map<std::string, std::string> data)
    : CollectionElementBase{std::move(data)}
{
}

std::string TemplatePrepKit::PartNumber() const { return Get(Element::PartNumber); }

TemplatePrepKit& TemplatePrepKit::PartNumber(std::string s)
{
    Set(Element::PartNumber, s);
    return *this;
}

bool TemplatePrepKit::TemplatePrepKit::HasPartNumber() const { return Has(Element::PartNumber); }

std::string TemplatePrepKit::LeftAdaptorSequence() const
{
    return Get(Element::LeftAdaptorSequence);
}

TemplatePrepKit& TemplatePrepKit::LeftAdaptorSequence(std::string s)
{
    Set(Element::LeftAdaptorSequence, s);
    return *this;
}

bool TemplatePrepKit::HasLeftAdaptorSequence() const { return Has(Element::LeftAdaptorSequence); }

std::string TemplatePrepKit::LeftPrimerSequence() const { return Get(Element::LeftPrimerSequence); }

TemplatePrepKit& TemplatePrepKit::LeftPrimerSequence(std::string s)
{
    Set(Element::LeftPrimerSequence, s);
    return *this;
}

bool TemplatePrepKit::HasLeftPrimerSequence() const { return Has(Element::LeftPrimerSequence); }

std::string TemplatePrepKit::RightAdaptorSequence() const
{
    return Get(Element::RightAdaptorSequence);
}

TemplatePrepKit& TemplatePrepKit::RightAdaptorSequence(std::string s)
{
    Set(Element::RightAdaptorSequence, s);
    return *this;
}

bool TemplatePrepKit::HasRightAdaptorSequence() const { return Has(Element::RightAdaptorSequence); }

std::string TemplatePrepKit::RightPrimerSequence() const
{
    return Get(Element::RightPrimerSequence);
}

TemplatePrepKit& TemplatePrepKit::RightPrimerSequence(std::string s)
{
    Set(Element::RightPrimerSequence, s);
    return *this;
}

bool TemplatePrepKit::HasRightPrimerSequence() const { return Has(Element::RightPrimerSequence); }

// ----------------------
// CollectionMetadata
// ----------------------

CollectionMetadata::CollectionMetadata() = default;

CollectionMetadata::CollectionMetadata(
    std::string subreadSetName,
    boost::optional<PacBio::BAM::AutomationParameters> automationParameters,
    boost::optional<PacBio::BAM::BindingKit> bindingKit,
    boost::optional<PacBio::BAM::ControlKit> controlKit,
    boost::optional<PacBio::BAM::SequencingKitPlate> sequencingKit,
    boost::optional<PacBio::BAM::TemplatePrepKit> templatePrepKit)
    : subreadSetName_{std::move(subreadSetName)}
    , automationParameters_{std::move(automationParameters)}
    , bindingKit_{std::move(bindingKit)}
    , controlKit_{std::move(controlKit)}
    , sequencingKit_{std::move(sequencingKit)}
    , templatePrepKit_{std::move(templatePrepKit)}
{
}

const std::string& CollectionMetadata::SubreadSetName() const { return subreadSetName_; }

const PacBio::BAM::AutomationParameters& CollectionMetadata::AutomationParameters() const
{
    return automationParameters_.get();
}

PacBio::BAM::AutomationParameters& CollectionMetadata::AutomationParameters()
{
    return automationParameters_.get();
}

CollectionMetadata& CollectionMetadata::AutomationParamters(
    PacBio::BAM::AutomationParameters params)
{
    automationParameters_ = std::move(params);
    return *this;
}

bool CollectionMetadata::HasAutomationParameters() const
{
    return automationParameters_.is_initialized();
}

const PacBio::BAM::BindingKit& CollectionMetadata::BindingKit() const { return bindingKit_.get(); }

PacBio::BAM::BindingKit& CollectionMetadata::BindingKit() { return bindingKit_.get(); }

CollectionMetadata& CollectionMetadata::BindingKit(PacBio::BAM::BindingKit kit)
{
    bindingKit_ = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasBindingKit() const { return bindingKit_.is_initialized(); }

const PacBio::BAM::ControlKit& CollectionMetadata::ControlKit() const { return controlKit_.get(); }

PacBio::BAM::ControlKit& CollectionMetadata::ControlKit() { return controlKit_.get(); }

CollectionMetadata& CollectionMetadata::ControlKit(PacBio::BAM::ControlKit kit)
{
    controlKit_ = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasControlKit() const { return controlKit_.is_initialized(); }

const PacBio::BAM::SequencingKitPlate& CollectionMetadata::SequencingKitPlate() const
{
    return sequencingKit_.get();
}
PacBio::BAM::SequencingKitPlate& CollectionMetadata::SequencingKitPlate()
{
    return sequencingKit_.get();
}

CollectionMetadata& CollectionMetadata::SequencingKitPlate(PacBio::BAM::SequencingKitPlate kit)
{
    sequencingKit_ = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasSequencingKitPlate() const { return sequencingKit_.is_initialized(); }

const PacBio::BAM::TemplatePrepKit& CollectionMetadata::TemplatePrepKit() const
{
    return templatePrepKit_.get();
}

PacBio::BAM::TemplatePrepKit& CollectionMetadata::TemplatePrepKit()
{
    return templatePrepKit_.get();
}

CollectionMetadata& CollectionMetadata::TemplatePrepKit(PacBio::BAM::TemplatePrepKit kit)
{
    templatePrepKit_ = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasTemplatePrepKit() const { return templatePrepKit_.is_initialized(); }

// ----------------------
// RunMetadata
// ----------------------

CollectionMetadata RunMetadata::Collection(const std::string& metadataXmlFn)
{
    return RunMetadataParser::LoadCollection(metadataXmlFn);
}

CollectionMetadata RunMetadata::Collection(std::istream& in)
{
    return RunMetadataParser::LoadCollection(in);
}

std::map<std::string, CollectionMetadata> RunMetadata::Collections(
    const std::string& runMetadataXmlFn)
{
    return RunMetadataParser::LoadCollections(runMetadataXmlFn);
}
std::map<std::string, CollectionMetadata> RunMetadata::Collections(std::istream& in)
{
    return RunMetadataParser::LoadCollections(in);
}

}  // namespace BAM
}  // namespace PacBio
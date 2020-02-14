// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/CollectionMetadata.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <boost/algorithm/string/predicate.hpp>

#include "DataSetUtils.h"
#include "RunMetadataParser.h"

namespace PacBio {
namespace BAM {
namespace {

boost::optional<ControlKit::CustomSequence> UpdateControlKitCache(const ControlKit& kit)
{
    if (!kit.HasChild("CustomSequence")) return boost::none;

    const std::string& customSeq = kit.ChildText("CustomSequence");
    const auto lines = [](const std::string& input) {
        std::vector<std::string> result;
        size_t pos = 0;
        size_t found = input.find("\\n");
        while (found != std::string::npos) {
            result.push_back(input.substr(pos, found - pos));
            pos = found + 2;  // "\n"
            found = input.find("\\n", pos);
        }
        result.push_back(input.substr(pos));  // store last
        return result;
    }(customSeq);

    if (lines.size() != 6)
        throw std::runtime_error{"[pbbam] run metadata ERROR: malformatted CustomSequence node"};

    return ControlKit::CustomSequence{lines.at(1), lines.at(3), lines.at(5)};
}

void UpdateControlKit(const boost::optional<ControlKit::CustomSequence>& cache, ControlKit& kit)
{
    std::ostringstream seq;
    seq << ">left_adapter\\n"
        << cache->LeftAdapter << "\\n"
        << ">right_adapter\\n"
        << cache->RightAdapter << "\\n"
        << ">custom_sequence\\n"
        << cache->Sequence;
    kit.ChildText("CustomSequence", seq.str());
}

}  // namespace

// ----------------------
// Automation
// ----------------------

Automation::Automation() : internal::DataSetElement{"Automation", XsdType::COLLECTION_METADATA} {}

Automation::Automation(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
}

DEFINE_ACCESSORS(Automation, AutomationParameters, AutomationParameters)

Automation& Automation::AutomationParameters(PacBio::BAM::AutomationParameters params)
{
    AutomationParameters() = params;
    return *this;
}
bool Automation::HasAutomationParameters() const { return HasChild(Element::AutomationParameters); }

// ----------------------
// AutomationParameter
// ----------------------

AutomationParameter::AutomationParameter()
    : internal::DataSetElement{"AutomationParameter", XsdType::BASE_DATA_MODEL}
{
}

AutomationParameter::AutomationParameter(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::BASE_DATA_MODEL}
{
}

AutomationParameter::AutomationParameter(const std::string& name, const std::string& type,
                                         const std::string& value)
    : internal::DataSetElement{"AutomationParameter", XsdType::BASE_DATA_MODEL}
{
    Name(name);
    Type(type);
    Value(value);
}

AutomationParameter::AutomationParameter(const std::string& name, const std::string& type,
                                         const std::string& value,
                                         const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
    Name(name);
    Type(type);
    Value(value);
}

const std::string& AutomationParameter::Name() const { return Attribute("Name"); }
std::string& AutomationParameter::Name() { return Attribute("Name"); }
AutomationParameter& AutomationParameter::Name(const std::string& name)
{
    Attribute("Name") = name;
    return *this;
}

const std::string& AutomationParameter::Type() const { return Attribute("ValueDataType"); }
std::string& AutomationParameter::Type() { return Attribute("ValueDataType"); }
AutomationParameter& AutomationParameter::Type(const std::string& type)
{
    Attribute("ValueDataType") = type;
    return *this;
}

const std::string& AutomationParameter::Value() const { return Attribute("SimpleValue"); }
std::string& AutomationParameter::Value() { return Attribute("SimpleValue"); }
AutomationParameter& AutomationParameter::Value(const std::string& value)
{
    Attribute("SimpleValue") = value;
    return *this;
}

// ----------------------
// AutomationParameters
// ----------------------

AutomationParameters::AutomationParameters()
    : internal::DataSetElement{"AutomationParameters", XsdType::BASE_DATA_MODEL}
{
}
AutomationParameters::AutomationParameters(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::BASE_DATA_MODEL}
{
}

AutomationParameters::iterator_type AutomationParameters::begin()
{
    return AutomationParameters::iterator_type(this, 0);
}

AutomationParameters::const_iterator_type AutomationParameters::begin() const { return cbegin(); }

AutomationParameters::const_iterator_type AutomationParameters::cbegin() const
{
    return AutomationParameters::const_iterator_type(this, 0);
}

AutomationParameters::iterator_type AutomationParameters::end()
{
    return AutomationParameters::iterator_type(this, NumChildren());
}

AutomationParameters::const_iterator_type AutomationParameters::end() const { return cend(); }

AutomationParameters::const_iterator_type AutomationParameters::cend() const
{
    return AutomationParameters::const_iterator_type(this, NumChildren());
}

int32_t AutomationParameters::CellNFCIndex() const
{
    return std::stoi(GetParameter(Element::CellNFCIndex));
}

AutomationParameters& AutomationParameters::CellNFCIndex(int32_t i)
{
    return SetParameter(Element::CellNFCIndex, "Int32", std::to_string(i));
}

bool AutomationParameters::HasCellNFCIndex() const { return HasParameter(Element::CellNFCIndex); }

int32_t AutomationParameters::CollectionNumber() const
{
    return std::stoi(GetParameter(Element::CollectionNumber));
}

AutomationParameters& AutomationParameters::CollectionNumber(int32_t i)
{
    return SetParameter(Element::CollectionNumber, "Int32", std::to_string(i));
}

bool AutomationParameters::HasCollectionNumber() const
{
    return HasParameter(Element::CellNFCIndex);
}

double AutomationParameters::Exposure() const { return std::stod(GetParameter(Element::Exposure)); }

AutomationParameters& AutomationParameters::Exposure(double d)
{
    return SetParameter(Element::Exposure, "Double", std::to_string(d));
}

bool AutomationParameters::HasExposure() const { return HasParameter(Element::Exposure); }

bool AutomationParameters::ExtendFirst() const
{
    return boost::iequals(GetParameter(Element::ExtendFirst), "True");
}

AutomationParameters& AutomationParameters::ExtendFirst(bool ok)
{
    return SetParameter(Element::ExtendFirst, "Boolean", (ok ? "True" : "False"));
}

bool AutomationParameters::HasExtendFirst() const { return HasParameter(Element::ExtendFirst); }

double AutomationParameters::ExtensionTime() const
{
    return std::stod(GetParameter(Element::ExtensionTime));
}

AutomationParameters& AutomationParameters::ExtensionTime(double d)
{
    return SetParameter(Element::ExtensionTime, "Double", std::to_string(d));
}

bool AutomationParameters::HasExtensionTime() const { return HasParameter(Element::ExtensionTime); }

int32_t AutomationParameters::ExtraIMWashes() const
{
    return std::stoi(GetParameter(Element::ExtraIMWashes));
}

AutomationParameters& AutomationParameters::ExtraIMWashes(int32_t i)
{
    return SetParameter(Element::ExtraIMWashes, "Int32", std::to_string(i));
}

bool AutomationParameters::HasExtraIMWashes() const { return HasParameter(Element::ExtraIMWashes); }

bool AutomationParameters::HasN2Switch() const
{
    return boost::iequals(GetParameter(Element::HasN2Switch), "True");
}

AutomationParameters& AutomationParameters::HasN2Switch(bool ok)
{
    return SetParameter(Element::HasN2Switch, "Boolean", (ok ? "True" : "False"));
}

bool AutomationParameters::HasHasN2Switch() const { return HasParameter(Element::HasN2Switch); }

std::string AutomationParameters::HQRFMethod() const { return GetParameter(Element::HQRFMethod); }

AutomationParameters& AutomationParameters::HQRFMethod(std::string s)
{
    return SetParameter(Element::HQRFMethod, "String", s);
}

bool AutomationParameters::HasHQRFMethod() const { return HasParameter(Element::HQRFMethod); }

double AutomationParameters::ImmobilizationTime() const
{
    return std::stod(GetParameter(Element::ImmobilizationTime));
}

AutomationParameters& AutomationParameters::ImmobilizationTime(double d)
{
    return SetParameter(Element::ImmobilizationTime, "Double", std::to_string(d));
}

bool AutomationParameters::HasImmobilizationTime() const
{
    return HasParameter(Element::ImmobilizationTime);
}

int32_t AutomationParameters::InsertSize() const
{
    return std::stoi(GetParameter(Element::InsertSize));
}

AutomationParameters& AutomationParameters::InsertSize(int32_t i)
{
    return SetParameter(Element::InsertSize, "Int32", std::to_string(i));
}

bool AutomationParameters::HasInsertSize() const { return HasParameter(Element::InsertSize); }

double AutomationParameters::MovieLength() const
{
    return std::stod(GetParameter(Element::MovieLength));
}

AutomationParameters& AutomationParameters::MovieLength(double d)
{
    return SetParameter(Element::ImmobilizationTime, "Double", std::to_string(d));
}

bool AutomationParameters::HasMovieLength() const { return HasParameter(Element::MovieLength); }

bool AutomationParameters::PCDinPlate() const
{
    return boost::iequals(GetParameter(Element::PCDinPlate), "True");
}

AutomationParameters& AutomationParameters::PCDinPlate(bool ok)
{
    return SetParameter(Element::PCDinPlate, "Boolean", (ok ? "True" : "False"));
}

bool AutomationParameters::HasPCDinPlate() const { return HasParameter(Element::PCDinPlate); }

bool AutomationParameters::PreExtensionWorkflow() const
{
    return boost::iequals(GetParameter(Element::PreExtensionWorkflow), "True");
}

AutomationParameters& AutomationParameters::PreExtensionWorkflow(bool ok)
{
    return SetParameter(Element::PreExtensionWorkflow, "Boolean", (ok ? "True" : "False"));
}

bool AutomationParameters::HasPreExtensionWorkflow() const
{
    return HasParameter(Element::PreExtensionWorkflow);
}

double AutomationParameters::SNRCut() const { return std::stod(GetParameter(Element::SNRCut)); }

AutomationParameters& AutomationParameters::SNRCut(double d)
{
    return SetParameter(Element::SNRCut, "Double", std::to_string(d));
}

bool AutomationParameters::HasSNRCut() const { return HasParameter(Element::SNRCut); }

int32_t AutomationParameters::TipSearchMaxDuration() const
{
    return std::stoi(GetParameter(Element::TipSearchMaxDuration));
}

AutomationParameters& AutomationParameters::TipSearchMaxDuration(int32_t i)
{
    return SetParameter(Element::TipSearchMaxDuration, "Int32", std::to_string(i));
}

bool AutomationParameters::HasTipSearchMaxDuration() const
{
    return HasParameter(Element::TipSearchMaxDuration);
}

bool AutomationParameters::UseStageHotStart() const
{
    return boost::iequals(GetParameter(Element::UseStageHotStart), "True");
}

AutomationParameters& AutomationParameters::UseStageHotStart(bool ok)
{
    return SetParameter(Element::UseStageHotStart, "Boolean", (ok ? "True" : "False"));
    return *this;
}

bool AutomationParameters::HasUseStageHotStart() const
{
    return HasParameter(Element::UseStageHotStart);
}

std::string AutomationParameters::GetParameter(const std::string& param) const
{
    const size_t count = NumChildren();
    for (size_t i = 0; i < count; ++i) {
        const internal::DataSetElement& child = *(children_.at(i).get());
        if (child.Attribute("Name") == param) return child.Attribute("SimpleValue");
    }

    throw std::runtime_error{""};
}

AutomationParameters& AutomationParameters::SetParameter(const std::string& name,
                                                         const std::string& type,
                                                         const std::string& value)
{
    const size_t count = NumChildren();
    for (size_t i = 0; i < count; ++i) {
        internal::DataSetElement* child = children_.at(i).get();
        if (child->Attribute("Name") == name) {
            child->Attribute("ValueDataType", type);
            child->Attribute("SimpleValue", value);
            return *this;
        }
    }

    // not found
    AddChild(AutomationParameter{name, type, value, internal::FromInputXml{}});
    return *this;
}

bool AutomationParameters::HasParameter(const std::string& param) const
{
    const size_t count = NumChildren();
    for (size_t i = 0; i < count; ++i) {
        const internal::DataSetElement* child = children_.at(i).get();
        if (child->Attribute("Name") == param) return true;
    }

    return false;
}

// ----------------------
// BindingKit
// ----------------------

BindingKit::BindingKit() : internal::DataSetElement{"BindingKit", XsdType::COLLECTION_METADATA} {}

BindingKit::BindingKit(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
}

const std::string& BindingKit::PartNumber() const { return Attribute(Element::PartNumber); }

std::string& BindingKit::PartNumber() { return Attribute(Element::PartNumber); }

BindingKit& BindingKit::PartNumber(std::string s)
{
    Attribute(Element::PartNumber, s);
    return *this;
}

bool BindingKit::HasPartNumber() { return HasAttribute(Element::PartNumber); }

// ----------------------
// Collections
// ----------------------

Collections::Collections() : internal::DataSetElement{"Collections", XsdType::NONE}
{
    Attribute("xmlns", "http://pacificbiosciences.com/PacBioCollectionMetadata.xsd");
}

Collections::Collections(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::NONE}
{
    Attribute("xmlns", "http://pacificbiosciences.com/PacBioCollectionMetadata.xsd");
}

// ----------------------
// ControlKit
// ----------------------

ControlKit::ControlKit() : internal::DataSetElement{"ControlKit", XsdType::COLLECTION_METADATA} {}

ControlKit::ControlKit(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
}

const std::string& ControlKit::PartNumber() const { return Attribute(Element::PartNumber); }

std::string& ControlKit::PartNumber() { return Attribute(Element::PartNumber); }

ControlKit& ControlKit::PartNumber(std::string s)
{
    Attribute(Element::PartNumber, s);
    return *this;
}
bool ControlKit::HasPartNumber() const { return HasAttribute(Element::PartNumber); }

const std::string& ControlKit::LeftAdapter() const
{
    if (!cache_) cache_ = UpdateControlKitCache(*this);
    return cache_->LeftAdapter;
}

ControlKit& ControlKit::LeftAdapter(std::string s)
{
    if (!cache_) cache_ = UpdateControlKitCache(*this);
    cache_->LeftAdapter = s;
    UpdateControlKit(cache_, *this);
    return *this;
}

bool ControlKit::HasLeftAdapter() const { return !LeftAdapter().empty(); }

const std::string& ControlKit::RightAdapter() const
{
    if (!cache_) cache_ = UpdateControlKitCache(*this);
    return cache_->RightAdapter;
}

ControlKit& ControlKit::RightAdapter(std::string s)
{
    if (!cache_) cache_ = UpdateControlKitCache(*this);
    cache_->RightAdapter = s;
    UpdateControlKit(cache_, *this);
    return *this;
}

bool ControlKit::HasRightAdapter() const { return !RightAdapter().empty(); }

const std::string& ControlKit::Sequence() const
{
    if (!cache_) cache_ = UpdateControlKitCache(*this);
    return cache_->Sequence;
}

ControlKit& ControlKit::Sequence(std::string s)
{
    if (!cache_) cache_ = UpdateControlKitCache(*this);
    cache_->Sequence = s;
    UpdateControlKit(cache_, *this);
    return *this;
}

bool ControlKit::HasSequence() const { return !Sequence().empty(); }

// ----------------------
// PPAConfig
// ----------------------

PPAConfig::PPAConfig() : internal::DataSetElement{"PPAConfig", XsdType::COLLECTION_METADATA} {}

PPAConfig::PPAConfig(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
}

const std::string& PPAConfig::Json() const { return Text(); }

std::string& PPAConfig::Json() { return Text(); }

PPAConfig& PPAConfig::Json(std::string json)
{
    Text(std::move(json));
    return *this;
}

// ----------------------
// SequencingKitPlate
// ----------------------

SequencingKitPlate::SequencingKitPlate()
    : internal::DataSetElement{"SequencingKitPlate", XsdType::COLLECTION_METADATA}
{
}

SequencingKitPlate::SequencingKitPlate(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
}

const std::string& SequencingKitPlate::PartNumber() const { return Attribute(Element::PartNumber); }

std::string& SequencingKitPlate::PartNumber() { return Attribute(Element::PartNumber); }

SequencingKitPlate& SequencingKitPlate::PartNumber(std::string s)
{
    Attribute(Element::PartNumber, s);
    return *this;
}

bool SequencingKitPlate::HasPartNumber() const { return HasAttribute(Element::PartNumber); }

// ----------------------
// TemplatePrepKit
// ----------------------

TemplatePrepKit::TemplatePrepKit()
    : internal::DataSetElement{"TemplatePrepKit", XsdType::COLLECTION_METADATA}
{
}

TemplatePrepKit::TemplatePrepKit(const internal::FromInputXml& fromInputXml)
    : internal::DataSetElement{"", fromInputXml, XsdType::COLLECTION_METADATA}
{
}

const std::string& TemplatePrepKit::PartNumber() const { return Attribute(Element::PartNumber); }

std::string& TemplatePrepKit::PartNumber() { return Attribute(Element::PartNumber); }

TemplatePrepKit& TemplatePrepKit::PartNumber(std::string s)
{
    Attribute(Element::PartNumber, s);
    return *this;
}

bool TemplatePrepKit::TemplatePrepKit::HasPartNumber() const
{
    return HasAttribute(Element::PartNumber);
}

std::string TemplatePrepKit::LeftAdaptorSequence() const
{
    return ChildText(Element::LeftAdaptorSequence);
}

TemplatePrepKit& TemplatePrepKit::LeftAdaptorSequence(std::string s)
{
    ChildText(Element::LeftAdaptorSequence, s);
    return *this;
}

bool TemplatePrepKit::HasLeftAdaptorSequence() const
{
    return HasChild(Element::LeftAdaptorSequence);
}

std::string TemplatePrepKit::LeftPrimerSequence() const
{
    return ChildText(Element::LeftPrimerSequence);
}

TemplatePrepKit& TemplatePrepKit::LeftPrimerSequence(std::string s)
{
    ChildText(Element::LeftPrimerSequence, s);
    return *this;
}

bool TemplatePrepKit::HasLeftPrimerSequence() const
{
    return HasChild(Element::LeftPrimerSequence);
}

std::string TemplatePrepKit::RightAdaptorSequence() const
{
    return ChildText(Element::RightAdaptorSequence);
}

TemplatePrepKit& TemplatePrepKit::RightAdaptorSequence(std::string s)
{
    ChildText(Element::RightAdaptorSequence, s);
    return *this;
}

bool TemplatePrepKit::HasRightAdaptorSequence() const
{
    return HasChild(Element::RightAdaptorSequence);
}

std::string TemplatePrepKit::RightPrimerSequence() const
{
    return ChildText(Element::RightPrimerSequence);
}

TemplatePrepKit& TemplatePrepKit::RightPrimerSequence(std::string s)
{
    ChildText(Element::RightPrimerSequence, s);
    return *this;
}

bool TemplatePrepKit::HasRightPrimerSequence() const
{
    return HasChild(Element::RightPrimerSequence);
}

// ----------------------
// CollectionMetadata
// ----------------------

CollectionMetadata::CollectionMetadata()
    : internal::StrictEntityType{"CollectionMetadata", "CollectionMetadata",
                                 XsdType::COLLECTION_METADATA}
{
}

CollectionMetadata::CollectionMetadata(const internal::FromInputXml& fromInputXml)
    : internal::StrictEntityType{"CollectionMetadata", "CollectionMetadata", fromInputXml,
                                 XsdType::COLLECTION_METADATA}
{
}

CollectionMetadata::CollectionMetadata(std::string subreadSetName)
    : internal::StrictEntityType{"CollectionMetadata", "CollectionMetadata",
                                 XsdType::COLLECTION_METADATA}
    , subreadSetName_{std::move(subreadSetName)}
{
}

CollectionMetadata::CollectionMetadata(std::string subreadSetName,
                                       const internal::FromInputXml& fromInputXml)
    : internal::StrictEntityType{"CollectionMetadata", "CollectionMetadata", fromInputXml,
                                 XsdType::COLLECTION_METADATA}
    , subreadSetName_{std::move(subreadSetName)}
{
}

const std::string& CollectionMetadata::SubreadSetName() const { return subreadSetName_; }

DEFINE_ACCESSORS(CollectionMetadata, Automation, Automation)

CollectionMetadata& CollectionMetadata::Automation(PacBio::BAM::Automation automation)
{
    Automation() = automation;
    return *this;
}

bool CollectionMetadata::HasAutomation() const { return HasChild(Element::Automation); }

const PacBio::BAM::AutomationParameters& CollectionMetadata::AutomationParameters() const
{
    const PacBio::BAM::Automation& automation = Automation();
    return automation.AutomationParameters();
}

PacBio::BAM::AutomationParameters& CollectionMetadata::AutomationParameters()
{
    PacBio::BAM::Automation& automation = Automation();
    return automation.AutomationParameters();
}

CollectionMetadata& CollectionMetadata::AutomationParameters(
    PacBio::BAM::AutomationParameters params)
{
    // PacBio::BAM::Automation& automation = Automation();
    AutomationParameters() = params;
    return *this;
}

bool CollectionMetadata::HasAutomationParameters() const
{
    return HasAutomation() && Automation().HasAutomationParameters();
}

DEFINE_ACCESSORS(CollectionMetadata, BindingKit, BindingKit)

CollectionMetadata& CollectionMetadata::BindingKit(PacBio::BAM::BindingKit kit)
{
    BindingKit() = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasBindingKit() const { return HasChild("BindingKit"); }

DEFINE_ACCESSORS(CollectionMetadata, ControlKit, ControlKit)

CollectionMetadata& CollectionMetadata::ControlKit(PacBio::BAM::ControlKit kit)
{
    ControlKit() = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasControlKit() const { return HasChild("ControlKit"); }

DEFINE_ACCESSORS(CollectionMetadata, PPAConfig, PPAConfig)

CollectionMetadata& CollectionMetadata::PPAConfig(PacBio::BAM::PPAConfig config)
{
    PPAConfig() = std::move(config);
    return *this;
}

bool CollectionMetadata::HasPPAConfig() const { return HasChild("PPAConfig"); }

DEFINE_ACCESSORS(CollectionMetadata, SequencingKitPlate, SequencingKitPlate)

CollectionMetadata& CollectionMetadata::SequencingKitPlate(PacBio::BAM::SequencingKitPlate kit)
{
    SequencingKitPlate() = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasSequencingKitPlate() const { return HasChild("SequencingKitPlate"); }

DEFINE_ACCESSORS(CollectionMetadata, TemplatePrepKit, TemplatePrepKit)

CollectionMetadata& CollectionMetadata::TemplatePrepKit(PacBio::BAM::TemplatePrepKit kit)
{
    TemplatePrepKit() = std::move(kit);
    return *this;
}

bool CollectionMetadata::HasTemplatePrepKit() const { return HasChild("TemplatePrepKit"); }

}  // namespace BAM
}  // namespace PacBio
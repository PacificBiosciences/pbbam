#ifndef PBBAM_COLLECTIONMETADATA_H
#define PBBAM_COLLECTIONMETADATA_H

#include <pbbam/Config.h>

#include <map>
#include <string>

#include <boost/optional.hpp>

#include <pbbam/internal/DataSetBaseTypes.h>

namespace PacBio {
namespace BAM {

class AutomationParameter : public internal::DataSetElement
{
public:
    AutomationParameter();
    AutomationParameter(const internal::FromInputXml& fromInputXml);
    AutomationParameter(const std::string& name, const std::string& type, const std::string& value);
    AutomationParameter(const std::string& name, const std::string& type, const std::string& value,
                        const internal::FromInputXml& fromInputXml);

    const std::string& Name() const;
    std::string& Name();
    AutomationParameter& Name(const std::string& name);

    const std::string& Type() const;
    std::string& Type();
    AutomationParameter& Type(const std::string& type);

    const std::string& Value() const;
    std::string& Value();
    AutomationParameter& Value(const std::string& value);
};

class AutomationParameters : public internal::DataSetElement
{
public:
    AutomationParameters();
    AutomationParameters(const internal::FromInputXml& fromInputXml);

public:
    int32_t CellNFCIndex() const;
    AutomationParameters& CellNFCIndex(int32_t i);
    bool HasCellNFCIndex() const;

    int32_t CollectionNumber() const;
    AutomationParameters& CollectionNumber(int32_t i);
    bool HasCollectionNumber() const;

    double Exposure() const;
    AutomationParameters& Exposure(double d);
    bool HasExposure() const;

    bool ExtendFirst() const;
    AutomationParameters& ExtendFirst(bool ok);
    bool HasExtendFirst() const;

    double ExtensionTime() const;
    AutomationParameters& ExtensionTime(double d);
    bool HasExtensionTime() const;

    int32_t ExtraIMWashes() const;
    AutomationParameters& ExtraIMWashes(int32_t i);
    bool HasExtraIMWashes() const;

    bool HasN2Switch() const;
    AutomationParameters& HasN2Switch(bool ok);
    bool HasHasN2Switch() const;

    std::string HQRFMethod() const;
    AutomationParameters& HQRFMethod(std::string s);
    bool HasHQRFMethod() const;

    double ImmobilizationTime() const;
    AutomationParameters& ImmobilizationTime(double d);
    bool HasImmobilizationTime() const;

    int32_t InsertSize() const;
    AutomationParameters& InsertSize(int32_t i);
    bool HasInsertSize() const;

    double MovieLength() const;
    AutomationParameters& MovieLength(double d);
    bool HasMovieLength() const;

    bool PCDinPlate() const;
    AutomationParameters& PCDinPlate(bool ok);
    bool HasPCDinPlate() const;

    bool PreExtensionWorkflow() const;
    AutomationParameters& PreExtensionWorkflow(bool ok);
    bool HasPreExtensionWorkflow() const;

    double SNRCut() const;
    AutomationParameters& SNRCut(double d);
    bool HasSNRCut() const;

    int32_t TipSearchMaxDuration() const;
    AutomationParameters& TipSearchMaxDuration(int32_t i);
    bool HasTipSearchMaxDuration() const;

    bool UseStageHotStart() const;
    AutomationParameters& UseStageHotStart(bool ok);
    bool HasUseStageHotStart() const;

public:
    // generic access & STL-iteration compatibility
    std::string GetParameter(const std::string& param) const;
    AutomationParameters& SetParameter(const std::string& name, const std::string& type,
                                       const std::string& value);
    bool HasParameter(const std::string& param) const;

    using value_type = AutomationParameter;
    using iterator_type = internal::DataSetElementIterator<value_type>;
    using const_iterator_type = internal::DataSetElementConstIterator<value_type>;

    const value_type& operator[](size_t index) const;
    value_type& operator[](size_t index);

    iterator_type begin();
    const_iterator_type begin() const;
    const_iterator_type cbegin() const;
    iterator_type end();
    const_iterator_type end() const;
    const_iterator_type cend() const;
};

class Automation : public internal::DataSetElement
{
public:
    Automation();
    Automation(const internal::FromInputXml& fromInputXml);

    const BAM::AutomationParameters& AutomationParameters() const;
    BAM::AutomationParameters& AutomationParameters();
    BAM::Automation& AutomationParameters(BAM::AutomationParameters params);
    bool HasAutomationParameters() const;
};

class BindingKit : public internal::DataSetElement
{
public:
    BindingKit();
    BindingKit(const internal::FromInputXml& fromInputXml);

    const std::string& PartNumber() const;
    std::string& PartNumber();
    BindingKit& PartNumber(std::string s);
    bool HasPartNumber();
};

class ControlKit : public internal::DataSetElement
{
public:
    ControlKit();
    ControlKit(const internal::FromInputXml& fromInputXml);

    const std::string& PartNumber() const;
    std::string& PartNumber();
    ControlKit& PartNumber(std::string s);
    bool HasPartNumber() const;

    const std::string& LeftAdapter() const;
    ControlKit& LeftAdapter(std::string s);
    bool HasLeftAdapter() const;

    const std::string& RightAdapter() const;
    ControlKit& RightAdapter(std::string s);
    bool HasRightAdapter() const;

    const std::string& Sequence() const;
    ControlKit& Sequence(std::string s);
    bool HasSequence() const;

    struct CustomSequence
    {
        std::string LeftAdapter;
        std::string RightAdapter;
        std::string Sequence;
    };

private:
    mutable boost::optional<CustomSequence> cache_;
};

class PPAConfig : public internal::DataSetElement
{
public:
    PPAConfig();
    PPAConfig(const internal::FromInputXml& fromInputXml);

    const std::string& Json() const;
    std::string& Json();
    PPAConfig& Json(std::string json);
};

class SequencingKitPlate : public internal::DataSetElement
{
public:
    SequencingKitPlate();
    SequencingKitPlate(const internal::FromInputXml& fromInputXml);

    const std::string& PartNumber() const;
    std::string& PartNumber();
    SequencingKitPlate& PartNumber(std::string s);
    bool HasPartNumber() const;
};

class TemplatePrepKit : public internal::DataSetElement
{
public:
    TemplatePrepKit();
    TemplatePrepKit(const internal::FromInputXml& fromInputXml);

    const std::string& PartNumber() const;
    std::string& PartNumber();
    TemplatePrepKit& PartNumber(std::string s);
    bool HasPartNumber() const;

    std::string LeftAdaptorSequence() const;
    TemplatePrepKit& LeftAdaptorSequence(std::string s);
    bool HasLeftAdaptorSequence() const;

    std::string LeftPrimerSequence() const;
    TemplatePrepKit& LeftPrimerSequence(std::string s);
    bool HasLeftPrimerSequence() const;

    std::string RightAdaptorSequence() const;
    TemplatePrepKit& RightAdaptorSequence(std::string s);
    bool HasRightAdaptorSequence() const;

    std::string RightPrimerSequence() const;
    TemplatePrepKit& RightPrimerSequence(std::string s);
    bool HasRightPrimerSequence() const;
};

class CollectionMetadata : public internal::StrictEntityType
{
public:
    static CollectionMetadata FromRawXml(const std::string& xml);

    CollectionMetadata();
    CollectionMetadata(const internal::FromInputXml& fromInputXml);
    CollectionMetadata(std::string subreadSetName);
    CollectionMetadata(std::string subreadSetName, const internal::FromInputXml& fromInputXml);

    const std::string& SubreadSetName() const;

    const BAM::Automation& Automation() const;
    BAM::Automation& Automation();
    CollectionMetadata& Automation(BAM::Automation automation);
    bool HasAutomation() const;

    const BAM::AutomationParameters& AutomationParameters() const;
    BAM::AutomationParameters& AutomationParameters();
    CollectionMetadata& AutomationParameters(BAM::AutomationParameters params);
    bool HasAutomationParameters() const;

    const BAM::BindingKit& BindingKit() const;
    BAM::BindingKit& BindingKit();
    CollectionMetadata& BindingKit(BAM::BindingKit kit);
    bool HasBindingKit() const;

    const BAM::ControlKit& ControlKit() const;
    BAM::ControlKit& ControlKit();
    CollectionMetadata& ControlKit(BAM::ControlKit kit);
    bool HasControlKit() const;

    const BAM::PPAConfig& PPAConfig() const;
    BAM::PPAConfig& PPAConfig();
    CollectionMetadata& PPAConfig(BAM::PPAConfig config);
    bool HasPPAConfig() const;

    const BAM::SequencingKitPlate& SequencingKitPlate() const;
    BAM::SequencingKitPlate& SequencingKitPlate();
    CollectionMetadata& SequencingKitPlate(BAM::SequencingKitPlate kit);
    bool HasSequencingKitPlate() const;

    const BAM::TemplatePrepKit& TemplatePrepKit() const;
    BAM::TemplatePrepKit& TemplatePrepKit();
    CollectionMetadata& TemplatePrepKit(BAM::TemplatePrepKit kit);
    bool HasTemplatePrepKit() const;

private:
    std::string subreadSetName_;
};

class Collections : public internal::DataSetElement
{
public:
    Collections();
    Collections(const internal::FromInputXml& fromInputXml);
};

}  // namespace BAM
}  // namespace PacBio

#endif  //  PBBAM_COLLECTIONMETADATA_H

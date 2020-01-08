// Author: Derek Barnett

#ifndef RUNMETADATA_H
#define RUNMETADATA_H

#include "pbbam/Config.h"

#include <cstdint>

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "pbbam/internal/DataSetBaseTypes.h"

namespace PacBio {
namespace BAM {

class CollectionElementBase
{
protected:
    CollectionElementBase();
    explicit CollectionElementBase(std::map<std::string, std::string> data);

    const std::string& Get(const std::string& name) const;
    std::string& Get(const std::string& name);
    bool Has(const std::string& name) const;
    CollectionElementBase& Set(const std::string& name, const std::string& value);

    std::map<std::string, std::string> data_;
};

class AutomationParameters : private CollectionElementBase
{
public:
    AutomationParameters();
    explicit AutomationParameters(std::map<std::string, std::string> data);

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
    AutomationParameters& SetParameter(const std::string& name, const std::string& value);
    bool HasParameter(const std::string& param) const;

    using iterator = std::map<std::string, std::string>::iterator;
    using const_iterator = std::map<std::string, std::string>::const_iterator;

    const_iterator cbegin() const { return data_.cbegin(); }
    const_iterator begin() const { return data_.begin(); }
    iterator begin() { return data_.begin(); }

    const_iterator cend() const { return data_.cend(); }
    const_iterator end() const { return data_.end(); }
    iterator end() { return data_.end(); }

public:
    // generic access
};

class BindingKit : private CollectionElementBase
{
public:
    BindingKit();
    explicit BindingKit(std::map<std::string, std::string> data);

    std::string PartNumber() const;
    BindingKit& PartNumber(std::string s);
    bool HasPartNumber();
};

class ControlKit : private CollectionElementBase
{
public:
    ControlKit();
    explicit ControlKit(std::map<std::string, std::string> data);

    std::string PartNumber() const;
    ControlKit& PartNumber(std::string s);
    bool HasPartNumber() const;

    std::string LeftAdapter() const;
    ControlKit& LeftAdapter(std::string s);
    bool HasLeftAdapter() const;

    std::string RightAdapter() const;
    ControlKit& RightAdapter(std::string s);
    bool HasRightAdapter() const;

    std::string Sequence() const;
    ControlKit& Sequence(std::string s);
    bool HasSequence() const;
};

class SequencingKitPlate : private CollectionElementBase
{
public:
    SequencingKitPlate();
    explicit SequencingKitPlate(std::map<std::string, std::string> data);

    std::string PartNumber() const;
    SequencingKitPlate& PartNumber(std::string s);
    bool HasPartNumber() const;
};

class TemplatePrepKit : private CollectionElementBase
{
public:
    TemplatePrepKit();
    explicit TemplatePrepKit(std::map<std::string, std::string> data);

    std::string PartNumber() const;
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

class CollectionMetadata
{
public:
    CollectionMetadata();
    CollectionMetadata(
        std::string subreadSetName,
        boost::optional<PacBio::BAM::AutomationParameters> automationParameters = boost::none,
        boost::optional<PacBio::BAM::BindingKit> bindingKit = boost::none,
        boost::optional<PacBio::BAM::ControlKit> controlKit = boost::none,
        boost::optional<PacBio::BAM::SequencingKitPlate> sequencingKit = boost::none,
        boost::optional<PacBio::BAM::TemplatePrepKit> templatePrepKit = boost::none);

    const std::string& SubreadSetName() const;

    const PacBio::BAM::AutomationParameters& AutomationParameters() const;
    PacBio::BAM::AutomationParameters& AutomationParameters();
    CollectionMetadata& AutomationParamters(PacBio::BAM::AutomationParameters params);
    bool HasAutomationParameters() const;

    const PacBio::BAM::BindingKit& BindingKit() const;
    PacBio::BAM::BindingKit& BindingKit();
    CollectionMetadata& BindingKit(PacBio::BAM::BindingKit kit);
    bool HasBindingKit() const;

    const PacBio::BAM::ControlKit& ControlKit() const;
    PacBio::BAM::ControlKit& ControlKit();
    CollectionMetadata& ControlKit(PacBio::BAM::ControlKit kit);
    bool HasControlKit() const;

    const PacBio::BAM::SequencingKitPlate& SequencingKitPlate() const;
    PacBio::BAM::SequencingKitPlate& SequencingKitPlate();
    CollectionMetadata& SequencingKitPlate(PacBio::BAM::SequencingKitPlate kit);
    bool HasSequencingKitPlate() const;

    const PacBio::BAM::TemplatePrepKit& TemplatePrepKit() const;
    PacBio::BAM::TemplatePrepKit& TemplatePrepKit();
    CollectionMetadata& TemplatePrepKit(PacBio::BAM::TemplatePrepKit kit);
    bool HasTemplatePrepKit() const;

private:
    std::string subreadSetName_;
    boost::optional<PacBio::BAM::AutomationParameters> automationParameters_;
    boost::optional<PacBio::BAM::BindingKit> bindingKit_;
    boost::optional<PacBio::BAM::ControlKit> controlKit_;
    boost::optional<PacBio::BAM::SequencingKitPlate> sequencingKit_;
    boost::optional<PacBio::BAM::TemplatePrepKit> templatePrepKit_;
};

class RunMetadata
{
public:
    ///
    /// \brief Read a single CollectionMetadata from XML file (or input stream).
    ///
    /// Intended for use with '<id>.metadata.xml'.
    ///
    /// \throw std::runtime_error on failure to parse XML for CollectionMetadata,
    ///                           or if the number of SubreadSets found is != 1.
    ///
    static CollectionMetadata Collection(const std::string& metadataXmlFn);
    static CollectionMetadata Collection(std::istream& in);

    ///
    /// \brief Read multiple CollectionMetadata objects from XML file (or input stream).
    ///
    /// Intended for use with multi-collection 'id.run.metadata.xml', but will
    /// work for single collection input.
    ///
    /// \throw std::runtime_error on failure to parse XML for CollectionMetadata(s)
    ///
    /// \returns map {subread set name => CollectionMetadata}
    ///
    static std::map<std::string, CollectionMetadata> Collections(
        const std::string& runMetadataXmlFn);
    static std::map<std::string, CollectionMetadata> Collections(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // RUNMETADATA_H

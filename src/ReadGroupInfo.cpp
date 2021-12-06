#include "PbbamInternalConfig.h"

#include <pbbam/ReadGroupInfo.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <iomanip>
#include <ios>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>

#include <boost/algorithm/cxx14/equal.hpp>
#include <boost/algorithm/string.hpp>

#include <pbbam/MD5.h>
#include <pbbam/SamTagCodec.h>
#include <pbbam/StringUtilities.h>

#include "ChemistryTable.h"

namespace PacBio {
namespace BAM {
namespace {

static const std::string sam_ID{"ID"};
static const std::string sam_CN{"CN"};
static const std::string sam_DS{"DS"};
static const std::string sam_DT{"DT"};
static const std::string sam_FO{"FO"};
static const std::string sam_KS{"KS"};
static const std::string sam_LB{"LB"};
static const std::string sam_PG{"PG"};
static const std::string sam_PI{"PI"};
static const std::string sam_PL{"PL"};
static const std::string sam_PM{"PM"};
static const std::string sam_PU{"PU"};
static const std::string sam_SM{"SM"};
static const std::string sam_BC{"BC"};

static const std::string feature_DQ{"DeletionQV"};
static const std::string feature_DT{"DeletionTag"};
static const std::string feature_IQ{"InsertionQV"};
static const std::string feature_MQ{"MergeQV"};
static const std::string feature_SQ{"SubstitutionQV"};
static const std::string feature_ST{"SubstitutionTag"};
static const std::string feature_IP{"Ipd"};
static const std::string feature_PW{"PulseWidth"};
static const std::string feature_PM{"PkMid"};
static const std::string feature_PA{"PkMean"};
static const std::string feature_PI{"PkMid2"};
static const std::string feature_PS{"PkMean2"};
static const std::string feature_LT{"Label"};
static const std::string feature_PQ{"LabelQV"};
static const std::string feature_PT{"AltLabel"};
static const std::string feature_PV{"AltLabelQV"};
static const std::string feature_PG{"PulseMergeQV"};
static const std::string feature_PC{"PulseCall"};
static const std::string feature_PD{"PrePulseFrames"};
static const std::string feature_PX{"PulseCallWidth"};
static const std::string feature_SF{"StartFrame"};
static const std::string feature_PE{"PulseExclusion"};

static const std::string token_RT{"READTYPE"};
static const std::string token_BK{"BINDINGKIT"};
static const std::string token_SK{"SEQUENCINGKIT"};
static const std::string token_BV{"BASECALLERVERSION"};
static const std::string token_FR{"FRAMERATEHZ"};
static const std::string token_CT{"CONTROL"};

static const std::string token_BF{"BarcodeFile"};
static const std::string token_BH{"BarcodeHash"};
static const std::string token_BC{"BarcodeCount"};
static const std::string token_BM{"BarcodeMode"};
static const std::string token_BQ{"BarcodeQuality"};

static const std::string codec_RAW{"Frames"};
static const std::string codec_V1{"CodecV1"};

static const std::string barcodemode_NONE{"None"};
static const std::string barcodemode_SYM{"Symmetric"};
static const std::string barcodemode_ASYM{"Asymmetric"};
static const std::string barcodemode_TAIL{"Tailed"};

static const std::string barcodequal_NONE{"None"};
static const std::string barcodequal_SCORE{"Score"};
static const std::string barcodequal_PROB{"Probability"};

static const std::string platformModelType_ASTRO{"ASTRO"};
static const std::string platformModelType_RS{"RS"};
static const std::string platformModelType_SEQUEL{"SEQUEL"};
static const std::string platformModelType_SEQUELII{"SEQUELII"};

std::string BaseFeatureName(const BaseFeature& feature)
{
    static const std::unordered_map<BaseFeature, std::string> lookup{
        {BaseFeature::DELETION_QV, feature_DQ},
        {BaseFeature::DELETION_TAG, feature_DT},
        {BaseFeature::INSERTION_QV, feature_IQ},
        {BaseFeature::MERGE_QV, feature_MQ},
        {BaseFeature::SUBSTITUTION_QV, feature_SQ},
        {BaseFeature::SUBSTITUTION_TAG, feature_ST},
        {BaseFeature::IPD, feature_IP},
        {BaseFeature::PULSE_WIDTH, feature_PW},
        {BaseFeature::PKMID, feature_PM},
        {BaseFeature::PKMEAN, feature_PA},
        {BaseFeature::PKMID2, feature_PI},
        {BaseFeature::PKMEAN2, feature_PS},
        {BaseFeature::LABEL_QV, feature_PQ},
        {BaseFeature::ALT_LABEL, feature_PT},
        {BaseFeature::ALT_LABEL_QV, feature_PV},
        {BaseFeature::PULSE_MERGE_QV, feature_PG},
        {BaseFeature::PULSE_CALL, feature_PC},
        {BaseFeature::PRE_PULSE_FRAMES, feature_PD},
        {BaseFeature::PULSE_CALL_WIDTH, feature_PX},
        {BaseFeature::START_FRAME, feature_SF},
        {BaseFeature::PULSE_EXCLUSION, feature_PE}};

    const auto found = lookup.find(feature);
    if (found != lookup.cend()) {
        return found->second;
    }
    throw std::runtime_error{"[pbbam] read group ERROR: unrecognized base feature"};
}

std::string FrameCodecName(const Data::FrameCodec& codec, const Data::FrameEncoder& encoder)
{
    switch (codec) {
        case Data::FrameCodec::RAW:
            return codec_RAW;
        case Data::FrameCodec::V1:
            return codec_V1;
        case Data::FrameCodec::V2:
            return encoder.Name();
        default:
            throw std::runtime_error{"[pbbam] read group ERROR: unrecognized frame codec"};
    }
}

std::string BarcodeModeName(const BarcodeModeType& mode)
{
    static const std::unordered_map<BarcodeModeType, std::string> lookup{
        {BarcodeModeType::NONE, barcodemode_NONE},
        {BarcodeModeType::SYMMETRIC, barcodemode_SYM},
        {BarcodeModeType::ASYMMETRIC, barcodemode_ASYM},
        {BarcodeModeType::TAILED, barcodemode_TAIL}};

    const auto found = lookup.find(mode);
    if (found != lookup.cend()) {
        return found->second;
    }
    throw std::runtime_error{"[pbbam] read group ERROR: unrecognized barcode mode type"};
}

std::string BarcodeQualityName(const BarcodeQualityType& type)
{
    static const std::unordered_map<BarcodeQualityType, std::string> lookup{
        {BarcodeQualityType::NONE, barcodequal_NONE},
        {BarcodeQualityType::SCORE, barcodequal_SCORE},
        {BarcodeQualityType::PROBABILITY, barcodequal_PROB}};

    const auto found = lookup.find(type);
    if (found != lookup.cend()) {
        return found->second;
    }
    throw std::runtime_error{"[pbbam] read group ERROR: unrecognized barcode quality type"};
}

std::string PlatformModelName(const PlatformModelType& type)
{
    static const std::unordered_map<PlatformModelType, std::string> lookup{
        {PlatformModelType::ASTRO, platformModelType_ASTRO},
        {PlatformModelType::RS, platformModelType_RS},
        {PlatformModelType::SEQUEL, platformModelType_SEQUEL},
        {PlatformModelType::SEQUELII, platformModelType_SEQUELII}};

    const auto found = lookup.find(type);
    if (found != lookup.cend()) {
        return found->second;
    }
    throw std::runtime_error{"[pbbam] read group ERROR: unrecognized platform model type"};
}

static const std::map<std::string, BaseFeature> nameToFeature{
    {feature_DQ, BaseFeature::DELETION_QV},
    {feature_DT, BaseFeature::DELETION_TAG},
    {feature_IQ, BaseFeature::INSERTION_QV},
    {feature_MQ, BaseFeature::MERGE_QV},
    {feature_SQ, BaseFeature::SUBSTITUTION_QV},
    {feature_ST, BaseFeature::SUBSTITUTION_TAG},
    {feature_IP, BaseFeature::IPD},
    {feature_PW, BaseFeature::PULSE_WIDTH},
    {feature_PM, BaseFeature::PKMID},
    {feature_PA, BaseFeature::PKMEAN},
    {feature_PI, BaseFeature::PKMID2},
    {feature_PS, BaseFeature::PKMEAN2},
    {feature_PQ, BaseFeature::LABEL_QV},
    {feature_PT, BaseFeature::ALT_LABEL},
    {feature_PV, BaseFeature::ALT_LABEL_QV},
    {feature_PC, BaseFeature::PULSE_CALL},
    {feature_PG, BaseFeature::PULSE_MERGE_QV},
    {feature_PD, BaseFeature::PRE_PULSE_FRAMES},
    {feature_PX, BaseFeature::PULSE_CALL_WIDTH},
    {feature_SF, BaseFeature::START_FRAME},
    {feature_PE, BaseFeature::PULSE_EXCLUSION}};

static const std::map<std::string, Data::FrameCodec> nameToCodec{{codec_RAW, Data::FrameCodec::RAW},
                                                                 {codec_V1, Data::FrameCodec::V1}};

static const std::map<std::string, BarcodeModeType> nameToBarcodeMode{
    {barcodemode_NONE, BarcodeModeType::NONE},
    {barcodemode_SYM, BarcodeModeType::SYMMETRIC},
    {barcodemode_ASYM, BarcodeModeType::ASYMMETRIC},
    {barcodemode_TAIL, BarcodeModeType::TAILED}};

static const std::map<std::string, BarcodeQualityType> nameToBarcodeQuality{
    {barcodequal_NONE, BarcodeQualityType::NONE},
    {barcodequal_SCORE, BarcodeQualityType::SCORE},
    {barcodequal_PROB, BarcodeQualityType::PROBABILITY}};

static const std::map<std::string, PlatformModelType> nameToPlatformModel{
    {platformModelType_ASTRO, PlatformModelType::ASTRO},
    {platformModelType_RS, PlatformModelType::RS},
    {platformModelType_SEQUEL, PlatformModelType::SEQUEL},
    {platformModelType_SEQUELII, PlatformModelType::SEQUELII}};

bool IsLikelyBarcodeKey(const std::string& name) { return name.find("Barcode") == 0; }

bool IsBaseFeature(const std::string& name)
{
    return nameToFeature.find(name) != nameToFeature.cend();
}

BaseFeature BaseFeatureFromName(const std::string& name) { return nameToFeature.at(name); }

Data::FrameCodec FrameCodecFromName(const std::string& name)
{
    const auto foundCodec = nameToCodec.find(name);
    if (foundCodec != nameToCodec.cend()) {
        return foundCodec->second;
    } else if (name.find("CodecV2") == 0) {
        return Data::FrameCodec::V2;
    }

    throw std::runtime_error{"[pbbam] read group ERROR: unknown codec name '" + name + "'"};
}

Data::FrameEncoder FrameEncoderFromName(const std::string& name)
{
    if (name.find("CodecV2") == 0) {
        const auto codecParts = BAM::Split(name, '/');
        assert(codecParts.size() == 3);
        const int exponentBits = std::stoi(codecParts[1]);
        const int mantissaBits = std::stoi(codecParts[2]);
        return Data::V2FrameEncoder{exponentBits, mantissaBits};
    } else {
        return Data::V1FrameEncoder{};
    }  // default
}

BarcodeModeType BarcodeModeFromName(const std::string& name) { return nameToBarcodeMode.at(name); }

BarcodeQualityType BarcodeQualityFromName(const std::string& name)
{
    return nameToBarcodeQuality.at(name);
}

PlatformModelType PlatformModelFromName(std::string name) { return nameToPlatformModel.at(name); }

}  // namespace

ReadGroupInfo::ReadGroupInfo(std::string baseId, std::pair<uint16_t, uint16_t> barcodes)
{
    std::ostringstream id;
    id << baseId << '/' << std::to_string(barcodes.first) << "--"
       << std::to_string(barcodes.second);
    id_ = id.str();
    baseId_ = std::move(baseId);
    barcodes_ = std::move(barcodes);
}

ReadGroupInfo::ReadGroupInfo() : readType_{"UNKNOWN"} {}

ReadGroupInfo::ReadGroupInfo(std::string id) : readType_{"UNKNOWN"} { Id(std::move(id)); }

ReadGroupInfo::ReadGroupInfo(std::string movieName, std::string readType)
    : ReadGroupInfo{ReadGroupInfoConfig{
          std::move(movieName), std::move(readType), PlatformModelType::SEQUEL, {}}}
{
}

ReadGroupInfo::ReadGroupInfo(std::string movieName, std::string readType,
                             std::pair<uint16_t, uint16_t> barcodes)
    : ReadGroupInfo{ReadGroupInfoConfig{
          std::move(movieName), std::move(readType), PlatformModelType::SEQUEL, std::move(barcodes),
      }}
{
}

ReadGroupInfo::ReadGroupInfo(std::string movieName, std::string readType,
                             PlatformModelType platform)
    : ReadGroupInfo{ReadGroupInfoConfig{std::move(movieName), std::move(readType), platform, {}}}
{
}

ReadGroupInfo::ReadGroupInfo(std::string movieName, std::string readType,
                             PlatformModelType platform, std::pair<uint16_t, uint16_t> barcodes)
    : ReadGroupInfo{ReadGroupInfoConfig{std::move(movieName), std::move(readType), platform,
                                        std::move(barcodes)}}
{
}

ReadGroupInfo::ReadGroupInfo(ReadGroupInfoConfig config)
    : movieName_{config.MovieName}, readType_{config.ReadType}
{
    if (config.Barcodes) {
        Id(MakeReadGroupId(config.MovieName, config.ReadType, *config.Barcodes));
        barcodes_ = std::move(config.Barcodes);
    } else {
        Id(MakeReadGroupId(config.MovieName, config.ReadType));
    }
    if (config.Platform) {
        platformModel_ = *config.Platform;
    }
}

bool ReadGroupInfo::operator==(const ReadGroupInfo& other) const noexcept
{
    const auto lhsFields = std::tie(
        id_, sequencingCenter_, date_, flowOrder_, keySequence_, library_, programs_,
        platformModel_, predictedInsertSize_, movieName_, sample_, readType_, bindingKit_,
        sequencingKit_, basecallerVersion_, frameRateHz_, control_, ipdCodec_, pulseWidthCodec_,
        hasBarcodeData_, barcodeFile_, barcodeHash_, barcodeCount_, barcodeMode_, barcodeQuality_);

    const auto rhsFields = std::tie(
        other.id_, other.sequencingCenter_, other.date_, other.flowOrder_, other.keySequence_,
        other.library_, other.programs_, other.platformModel_, other.predictedInsertSize_,
        other.movieName_, other.sample_, other.readType_, other.bindingKit_, other.sequencingKit_,
        other.basecallerVersion_, other.frameRateHz_, other.control_, other.ipdCodec_,
        other.pulseWidthCodec_, other.hasBarcodeData_, other.barcodeFile_, other.barcodeHash_,
        other.barcodeCount_, other.barcodeMode_, other.barcodeQuality_);

    return lhsFields == rhsFields &&
           boost::algorithm::equal(features_.cbegin(), features_.cend(), other.features_.cbegin(),
                                   other.features_.cend()) &&
           boost::algorithm::equal(custom_.cbegin(), custom_.cend(), other.custom_.cbegin(),
                                   other.custom_.cend());
}

bool ReadGroupInfo::operator<(const ReadGroupInfo& other) const noexcept { return id_ < other.id_; }

size_t ReadGroupInfo::BarcodeCount() const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{
            "[pbbam] read group ERROR: barcode count requested but barcode data is missing"};
    }
    return barcodeCount_;
}

ReadGroupInfo& ReadGroupInfo::BarcodeData(std::string barcodeFile, std::string barcodeHash,
                                          size_t barcodeCount, BarcodeModeType barcodeMode,
                                          BarcodeQualityType barcodeQuality)
{
    barcodeFile_ = std::move(barcodeFile);
    barcodeHash_ = std::move(barcodeHash);
    barcodeCount_ = barcodeCount;
    barcodeMode_ = barcodeMode;
    barcodeQuality_ = barcodeQuality;
    hasBarcodeData_ = true;
    return *this;
}

std::string ReadGroupInfo::BarcodeFile() const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{"[pbbam] read group ERROR: barcode file field is missing"};
    }
    return barcodeFile_;
}

std::string ReadGroupInfo::BarcodeHash() const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{"[pbbam] read group ERROR: barcode hash field is missing"};
    }
    return barcodeHash_;
}

BarcodeModeType ReadGroupInfo::BarcodeMode() const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{"[pbbam] read group ERROR: barcode mode field is missing"};
    }
    return barcodeMode_;
}

BarcodeQualityType ReadGroupInfo::BarcodeQuality() const
{
    if (!hasBarcodeData_) {
        throw std::runtime_error{"[pbbam] read group ERROR: barcode quality field is missing"};
    }
    return barcodeQuality_;
}

boost::optional<uint16_t> ReadGroupInfo::BarcodeForward() const
{
    const auto barcodes = Barcodes();
    if (barcodes) {
        return barcodes->first;
    }
    return {};
}

boost::optional<uint16_t> ReadGroupInfo::BarcodeReverse() const
{
    const auto barcodes = Barcodes();
    if (barcodes) {
        return barcodes->second;
    }
    return {};
}

boost::optional<std::pair<uint16_t, uint16_t>> ReadGroupInfo::Barcodes() const { return barcodes_; }

std::string ReadGroupInfo::BarcodeSequence() const
{
    const auto found = custom_.find(sam_BC);
    if (found == custom_.cend()) {
        return {};
    }
    return found->second;
}

ReadGroupInfo& ReadGroupInfo::BarcodeSequence(std::string barcodeSequence)
{
    custom_[sam_BC] = std::move(barcodeSequence);
    return *this;
}

std::string ReadGroupInfo::BasecallerVersion() const { return basecallerVersion_; }

ReadGroupInfo& ReadGroupInfo::BasecallerVersion(std::string versionNumber)
{
    if (basecallerVersion_ != versionNumber) {
        basecallerVersion_ = std::move(versionNumber);
        sequencingChemistry_.clear();  // reset cached chemistry name
    }
    return *this;
}

std::string ReadGroupInfo::BaseFeatureTag(BaseFeature feature) const
{
    const auto iter = features_.find(feature);
    if (iter == features_.end()) {
        return {};
    }
    return iter->second;
}

ReadGroupInfo& ReadGroupInfo::BaseFeatureTag(BaseFeature feature, std::string tag)
{
    features_[feature] = std::move(tag);
    return *this;
}

std::string ReadGroupInfo::BaseId() const { return baseId_; }

std::string ReadGroupInfo::BindingKit() const { return bindingKit_; }

ReadGroupInfo& ReadGroupInfo::BindingKit(std::string kitNumber)
{
    if (bindingKit_ != kitNumber) {
        bindingKit_ = std::move(kitNumber);
        sequencingChemistry_.clear();  // reset cached chemistry name
    }
    return *this;
}

ReadGroupInfo& ReadGroupInfo::ClearBarcodeData()
{
    barcodeFile_.clear();
    barcodeHash_.clear();
    hasBarcodeData_ = false;
    return *this;
}

ReadGroupInfo& ReadGroupInfo::ClearBaseFeatures()
{
    features_.clear();
    return *this;
}

bool ReadGroupInfo::Control() const { return control_; }

ReadGroupInfo& ReadGroupInfo::Control(bool ctrl)
{
    control_ = ctrl;
    return *this;
}

std::map<std::string, std::string> ReadGroupInfo::CustomTags() const { return custom_; }

ReadGroupInfo& ReadGroupInfo::CustomTags(std::map<std::string, std::string> custom)
{
    custom_ = std::move(custom);
    return *this;
}

std::string ReadGroupInfo::Date() const { return date_; }

ReadGroupInfo& ReadGroupInfo::Date(std::string date)
{
    date_ = std::move(date);
    return *this;
}

void ReadGroupInfo::DecodeBarcodeKey(const std::string& key, std::string value)
{
    if (key == token_BF) {
        barcodeFile_ = std::move(value);
    } else if (key == token_BH) {
        barcodeHash_ = std::move(value);
    } else if (key == token_BC) {
        barcodeCount_ = std::stoul(value);
    } else if (key == token_BM) {
        barcodeMode_ = BarcodeModeFromName(value);
    } else if (key == token_BQ) {
        barcodeQuality_ = BarcodeQualityFromName(value);
    }
}

void ReadGroupInfo::DecodeFrameCodecKey(const std::string& key, std::string value)
{
    const auto keyParts = Split(key, ':');
    if (keyParts.size() == 2) {
        const auto& subkey = keyParts.at(0);
        if (subkey == feature_IP) {
            ipdCodec_ = FrameCodecFromName(keyParts.at(1));
            ipdEncoder_ = FrameEncoderFromName(keyParts.at(1));
            features_[BaseFeature::IPD] = std::move(value);
        } else if (subkey == feature_PW) {
            pulseWidthCodec_ = FrameCodecFromName(keyParts.at(1));
            pulseWidthEncoder_ = FrameEncoderFromName(keyParts.at(1));
            features_[BaseFeature::PULSE_WIDTH] = std::move(value);
        }
    }
}

void ReadGroupInfo::DecodeSamDescription(const std::string& description)
{
    const auto tokens = Split(description, ';');
    if (tokens.empty()) {
        return;
    }

    // iterate over tokens
    for (const auto& token : tokens) {

        const auto foundEqual = token.find('=');
        if (foundEqual == std::string::npos) {
            continue;
        }

        const auto key = token.substr(0, foundEqual);
        auto value = token.substr(foundEqual + 1);

        // 'mandatory' items
        if (key == token_RT) {
            readType_ = std::move(value);
        } else if (key == token_BK) {
            bindingKit_ = std::move(value);
        } else if (key == token_BV) {
            basecallerVersion_ = std::move(value);
        } else if (key == token_SK) {
            sequencingKit_ = std::move(value);
        } else if (key == token_FR) {
            frameRateHz_ = std::move(value);
        } else if (key == token_CT) {
            control_ = (value == "TRUE");

            // base features
        } else if (IsBaseFeature(key)) {
            features_[BaseFeatureFromName(key)] = std::move(value);

            // barcode data
        } else if (IsLikelyBarcodeKey(key)) {
            DecodeBarcodeKey(key, std::move(value));

            // frame codecs
        } else {
            DecodeFrameCodecKey(key, std::move(value));
        }
    }

    hasBarcodeData_ = !barcodeFile_.empty();
}

std::string ReadGroupInfo::EncodeSamDescription() const
{
    constexpr static const char SEP = ';';
    constexpr static const char COLON = ':';
    constexpr static const char EQ = '=';

    std::string result{token_RT + EQ + readType_};

    std::string featureName;
    for (const auto& feature : features_) {

        featureName = BaseFeatureName(feature.first);
        if (featureName.empty() || feature.second.empty()) {
            continue;
        } else if (featureName == feature_IP) {
            featureName.push_back(COLON);
            featureName.append(FrameCodecName(ipdCodec_, ipdEncoder_));
        } else if (featureName == feature_PW) {
            featureName.push_back(COLON);
            featureName.append(FrameCodecName(pulseWidthCodec_, pulseWidthEncoder_));
        }
        result.append(SEP + featureName + EQ + feature.second);
    }

    if (!bindingKit_.empty()) {
        result.append(SEP + token_BK + EQ + bindingKit_);
    }
    if (!sequencingKit_.empty()) {
        result.append(SEP + token_SK + EQ + sequencingKit_);
    }
    if (!basecallerVersion_.empty()) {
        result.append(SEP + token_BV + EQ + basecallerVersion_);
    }
    if (!frameRateHz_.empty()) {
        result.append(SEP + token_FR + EQ + frameRateHz_);
    }
    if (control_) {
        result.append(SEP + token_CT + EQ + (control_ ? "TRUE" : "FALSE"));
    }

    if (hasBarcodeData_) {
        const std::string barcodeData{SEP + token_BF + EQ + barcodeFile_ + SEP + token_BH + EQ +
                                      barcodeHash_ + SEP + token_BC + EQ +
                                      std::to_string(barcodeCount_) + SEP + token_BM + EQ +
                                      BarcodeModeName(barcodeMode_) + SEP + token_BQ + EQ +
                                      BarcodeQualityName(barcodeQuality_)};
        result.append(barcodeData);
    }

    return result;
}

std::string ReadGroupInfo::FlowOrder() const { return flowOrder_; }

ReadGroupInfo& ReadGroupInfo::FlowOrder(std::string order)
{
    flowOrder_ = std::move(order);
    return *this;
}

std::string ReadGroupInfo::FrameRateHz() const { return frameRateHz_; }

ReadGroupInfo& ReadGroupInfo::FrameRateHz(std::string frameRateHz)
{
    frameRateHz_ = std::move(frameRateHz);
    return *this;
}

ReadGroupInfo ReadGroupInfo::FromSam(const std::string& sam)
{
    // pop off '@RG\t', then split rest of line into tokens
    const auto tokens = Split(sam.substr(4), '\t');
    if (tokens.empty()) {
        return {};
    }

    ReadGroupInfo rg;
    std::map<std::string, std::string> custom;

    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set read group info
        if (tokenTag == sam_ID) {
            rg.Id(std::move(tokenValue));
        } else if (tokenTag == sam_CN) {
            rg.SequencingCenter(std::move(tokenValue));
        } else if (tokenTag == sam_DT) {
            rg.Date(std::move(tokenValue));
        } else if (tokenTag == sam_FO) {
            rg.FlowOrder(std::move(tokenValue));
        } else if (tokenTag == sam_KS) {
            rg.KeySequence(std::move(tokenValue));
        } else if (tokenTag == sam_LB) {
            rg.Library(std::move(tokenValue));
        } else if (tokenTag == sam_PG) {
            rg.Programs(std::move(tokenValue));
        } else if (tokenTag == sam_PI) {
            rg.PredictedInsertSize(std::move(tokenValue));
        } else if (tokenTag == sam_PU) {
            rg.MovieName(std::move(tokenValue));
        } else if (tokenTag == sam_SM) {
            rg.Sample(std::move(tokenValue));
        } else if (tokenTag == sam_DS) {
            rg.DecodeSamDescription(std::move(tokenValue));
        } else if (tokenTag == sam_PM) {
            rg.PlatformModel(PlatformModelFromName(std::move(tokenValue)));

            // if not platform name (always "PACBIO" for us), store as a custom tag
        } else if (tokenTag != sam_PL) {
            custom[tokenTag] = std::move(tokenValue);
        }
    }
    rg.CustomTags(std::move(custom));

    return rg;
}

std::string ReadGroupInfo::GetBaseId(const std::string& id)
{
    const auto slashAt = id.find('/');
    if (slashAt == std::string::npos) {
        return id;
    } else {
        return id.substr(0, slashAt);
    }
}

bool ReadGroupInfo::HasBarcodeData() const { return hasBarcodeData_; }

bool ReadGroupInfo::HasBaseFeature(BaseFeature feature) const
{
    return features_.find(feature) != features_.end();
}

std::string ReadGroupInfo::Id() const { return id_; }

ReadGroupInfo& ReadGroupInfo::Id(const std::string& movieName, const std::string& readType)
{
    return Id(MakeReadGroupId(movieName, readType));
}

ReadGroupInfo& ReadGroupInfo::Id(std::string id)
{
    barcodes_.reset();

    // maybe parse for barcode labels
    const auto slashAt = id.find('/');
    if (slashAt != std::string::npos) {
        // looks like we do, parse & store
        const auto tokens = Split(id.substr(slashAt + 1), '-');
        if (tokens.size() != 3) {
            throw std::runtime_error{
                "[pbbam] read group ERROR: could not fetch barcodes from malformed read group "
                "ID: " +
                id + " Must be in the form: {RGID_STRING}/{bcForward}--{bcReverse}"};
        }

        // catch here so we can give more informative message
        try {
            barcodes_ = std::pair<uint16_t, uint16_t>(static_cast<uint16_t>(std::stoul(tokens[0])),
                                                      static_cast<uint16_t>(std::stoul(tokens[2])));
        } catch (std::exception& e) {
            throw std::runtime_error{
                "[pbbam] read group ERROR: could not fetch barcodes from malformed read group "
                "ID: " +
                id + " Must be in the form: {RGID_STRING}/{bcForward}--{bcReverse}"};
        }
    }

    baseId_ = id.substr(0, slashAt);
    id_ = std::move(id);
    return *this;
}

int32_t ReadGroupInfo::IdToInt(const std::string& rgId)
{
    const auto id = GetBaseId(rgId);
    const uint32_t rawid = std::stoul(id, nullptr, 16);
    return static_cast<int32_t>(rawid);
}

std::string ReadGroupInfo::IntToId(const int32_t id)
{
    std::ostringstream s;
    s << std::setfill('0') << std::setw(8) << std::hex << id;
    return s.str();
}

Data::FrameCodec ReadGroupInfo::IpdCodec() const { return ipdCodec_; }

ReadGroupInfo& ReadGroupInfo::IpdCodec(Data::FrameCodec codec, std::string tag)
{
    // store desired codec type
    ipdCodec_ = std::move(codec);

    // update base features map
    const std::string actualTag = (tag.empty() ? "ip" : std::move(tag));
    BaseFeatureTag(BaseFeature::IPD, actualTag);
    return *this;
}

Data::FrameEncoder ReadGroupInfo::IpdFrameEncoder() const { return ipdEncoder_; }

ReadGroupInfo& ReadGroupInfo::IpdFrameEncoder(Data::FrameEncoder encoder)
{
    ipdEncoder_ = std::move(encoder);
    return *this;
}

bool ReadGroupInfo::IsValid() const { return !id_.empty(); }

std::string ReadGroupInfo::KeySequence() const { return keySequence_; }

ReadGroupInfo& ReadGroupInfo::KeySequence(std::string sequence)
{
    keySequence_ = std::move(sequence);
    return *this;
}

std::string ReadGroupInfo::Library() const { return library_; }

ReadGroupInfo& ReadGroupInfo::Library(std::string library)
{
    library_ = std::move(library);
    return *this;
}

std::string ReadGroupInfo::MovieName() const { return movieName_; }

ReadGroupInfo& ReadGroupInfo::MovieName(std::string movieName)
{
    movieName_ = std::move(movieName);
    return *this;
}

std::string ReadGroupInfo::Platform() const { return std::string("PACBIO"); }

PlatformModelType ReadGroupInfo::PlatformModel() const { return platformModel_; }

ReadGroupInfo& ReadGroupInfo::PlatformModel(PlatformModelType platform)
{
    platformModel_ = platform;
    return *this;
}

std::string ReadGroupInfo::PredictedInsertSize() const { return predictedInsertSize_; }

ReadGroupInfo& ReadGroupInfo::PredictedInsertSize(std::string size)
{
    predictedInsertSize_ = std::move(size);
    return *this;
}

std::string ReadGroupInfo::Programs() const { return programs_; }

ReadGroupInfo& ReadGroupInfo::Programs(std::string programs)
{
    programs_ = std::move(programs);
    return *this;
}

Data::FrameCodec ReadGroupInfo::PulseWidthCodec() const { return pulseWidthCodec_; }

ReadGroupInfo& ReadGroupInfo::PulseWidthCodec(Data::FrameCodec codec, std::string tag)
{
    // store desired codec type
    pulseWidthCodec_ = std::move(codec);

    // update base features map
    const std::string actualTag = (tag.empty() ? "pw" : std::move(tag));
    BaseFeatureTag(BaseFeature::PULSE_WIDTH, actualTag);
    return *this;
}

Data::FrameEncoder ReadGroupInfo::PulseWidthFrameEncoder() const { return pulseWidthEncoder_; }

ReadGroupInfo& ReadGroupInfo::PulseWidthFrameEncoder(Data::FrameEncoder encoder)
{
    pulseWidthEncoder_ = std::move(encoder);
    return *this;
}

std::string ReadGroupInfo::ReadType() const { return readType_; }

ReadGroupInfo& ReadGroupInfo::ReadType(std::string type)
{
    readType_ = std::move(type);
    return *this;
}

ReadGroupInfo& ReadGroupInfo::RemoveBaseFeature(BaseFeature feature)
{
    const auto iter = features_.find(feature);
    if (iter != features_.end()) {
        features_.erase(iter);
    }
    return *this;
}

std::string ReadGroupInfo::Sample() const { return sample_; }

ReadGroupInfo& ReadGroupInfo::Sample(std::string sample)
{
    sample_ = std::move(sample);
    return *this;
}

std::string ReadGroupInfo::SequencingCenter() const { return sequencingCenter_; }

ReadGroupInfo& ReadGroupInfo::SequencingCenter(std::string center)
{
    sequencingCenter_ = std::move(center);
    return *this;
}

std::string ReadGroupInfo::SequencingChemistry() const
{
    if (!sequencingChemistry_.empty()) {
        return sequencingChemistry_;
    }
    return sequencingChemistry_ =
               SequencingChemistryFromTriple(BindingKit(), SequencingKit(), BasecallerVersion());
}

std::string ReadGroupInfo::SequencingChemistryFromTriple(const std::string& bindingKit,
                                                         const std::string& sequencingKit,
                                                         const std::string& basecallerVersion)
{
    const auto verFields = Split(basecallerVersion, '.');
    if (verFields.size() < 2) {
        throw std::runtime_error{"[pbbam] read group ERROR: basecaller version is too short: " +
                                 basecallerVersion};
    }
    const std::string version{verFields.at(0) + '.' + verFields.at(1)};

    // check updated table first, if it exists (empty if not), overriding the built-in lookup
    for (const auto& row : GetChemistryTableFromEnv()) {
        if (bindingKit == row[0] && sequencingKit == row[1] && version == row[2]) {
            return row[3];
        }
    }

    for (const auto& row : BuiltInChemistryTable()) {
        if (bindingKit == row[0] && sequencingKit == row[1] && version == row[2]) {
            return row[3];
        }
    }

    // not found
    throw InvalidSequencingChemistryException{bindingKit, sequencingKit, basecallerVersion};
}

std::string ReadGroupInfo::SequencingKit() const { return sequencingKit_; }

ReadGroupInfo& ReadGroupInfo::SequencingKit(std::string kitNumber)
{
    if (sequencingKit_ != kitNumber) {
        sequencingKit_ = std::move(kitNumber);
        sequencingChemistry_.clear();  // reset cached chemistry name
    }
    return *this;
}

std::string ReadGroupInfo::ToSam(const ReadGroupInfo& rg) { return rg.ToSam(); }

std::string ReadGroupInfo::ToSam() const
{
    std::ostringstream out;
    out << "@RG" << MakeSamTag(sam_ID, id_) << MakeSamTag(sam_PL, Platform());

    const auto description = EncodeSamDescription();
    if (!description.empty()) {
        out << MakeSamTag(sam_DS, description);
    }

    if (!sequencingCenter_.empty()) {
        out << MakeSamTag(sam_CN, sequencingCenter_);
    }
    if (!date_.empty()) {
        out << MakeSamTag(sam_DT, date_);
    }
    if (!flowOrder_.empty()) {
        out << MakeSamTag(sam_FO, flowOrder_);
    }
    if (!keySequence_.empty()) {
        out << MakeSamTag(sam_KS, keySequence_);
    }
    if (!library_.empty()) {
        out << MakeSamTag(sam_LB, library_);
    }
    if (!programs_.empty()) {
        out << MakeSamTag(sam_PG, programs_);
    }
    if (!predictedInsertSize_.empty()) {
        out << MakeSamTag(sam_PI, predictedInsertSize_);
    }
    if (!movieName_.empty()) {
        out << MakeSamTag(sam_PU, movieName_);
    }
    if (!sample_.empty()) {
        out << MakeSamTag(sam_SM, sample_);
    }

    out << MakeSamTag(sam_PM, PlatformModelName(platformModel_));

    // append any custom tags
    for (const auto& attribute : custom_) {
        out << MakeSamTag(attribute.first, attribute.second);
    }

    return out.str();
}

// ---------------------------------------------------------

std::string MakeReadGroupId(const std::string& movieName, const std::string& readType)
{
    return MD5Hash(movieName + "//" + readType).substr(0, 8);
}

std::string MakeReadGroupId(const std::string& movieName, const std::string& readType,
                            const std::string& barcodeString)
{
    const std::string baseId{MakeReadGroupId(movieName, readType)};
    return baseId + "/" + barcodeString;
}

std::string MakeReadGroupId(const std::string& movieName, const std::string& readType,
                            const std::pair<int16_t, int16_t>& barcodes)
{
    const std::string barcodeString{std::to_string(barcodes.first) + "--" +
                                    std::to_string(barcodes.second)};
    return MakeReadGroupId(movieName, readType, barcodeString);
}

std::string MakeReadGroupId(const ReadGroupInfo& readGroup)
{
    const auto barcodes = readGroup.Barcodes();
    if (barcodes) {
        const int16_t bcFor = barcodes->first;
        const int16_t bcRev = barcodes->second;
        return MakeReadGroupId(readGroup.MovieName(), readGroup.ReadType(),
                               std::make_pair(bcFor, bcRev));
    } else {
        return MakeReadGroupId(readGroup.MovieName(), readGroup.ReadType());
    }
}

std::string MakeLegacyReadGroupId(const std::string& movieName, const std::string& readType)
{
    return MD5Hash(movieName + "//" + readType).substr(0, 8);
}

std::string MakeLegacyReadGroupId(const std::string& movieName, const std::string& readType,
                                  const std::string& barcodeString)
{
    const std::string baseId{
        MD5Hash(movieName + "//" + readType + "//" + barcodeString).substr(0, 8)};
    return baseId + "/" + barcodeString;
}

std::string MakeLegacyReadGroupId(const std::string& movieName, const std::string& readType,
                                  const std::pair<int16_t, int16_t>& barcodes)
{
    const std::string barcodeString{std::to_string(barcodes.first) + "--" +
                                    std::to_string(barcodes.second)};
    return MakeLegacyReadGroupId(movieName, readType, barcodeString);
}

std::string MakeLegacyReadGroupId(const ReadGroupInfo& readGroup)
{
    const auto barcodes = readGroup.Barcodes();
    if (barcodes) {
        const int16_t bcFor = barcodes->first;
        const int16_t bcRev = barcodes->second;
        return MakeLegacyReadGroupId(readGroup.MovieName(), readGroup.ReadType(),
                                     std::make_pair(bcFor, bcRev));
    } else {
        return MakeLegacyReadGroupId(readGroup.MovieName(), readGroup.ReadType());
    }
}

}  // namespace BAM
}  // namespace PacBio

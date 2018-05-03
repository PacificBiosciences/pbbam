// File Description
/// \file ReadGroupInfo.cpp
/// \brief Implements the ReadGroupInfo class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ReadGroupInfo.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

#include "ChemistryTable.h"
#include "SequenceUtils.h"
#include "pbbam/MD5.h"

namespace PacBio {
namespace BAM {
namespace internal {

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

// clang-format off
static std::string BaseFeatureName(const BaseFeature& feature)
{
    switch (feature) {
        case BaseFeature::DELETION_QV      : return feature_DQ;
        case BaseFeature::DELETION_TAG     : return feature_DT;
        case BaseFeature::INSERTION_QV     : return feature_IQ;
        case BaseFeature::MERGE_QV         : return feature_MQ;
        case BaseFeature::SUBSTITUTION_QV  : return feature_SQ;
        case BaseFeature::SUBSTITUTION_TAG : return feature_ST;
        case BaseFeature::IPD              : return feature_IP;
        case BaseFeature::PULSE_WIDTH      : return feature_PW;
        case BaseFeature::PKMID            : return feature_PM;
        case BaseFeature::PKMEAN           : return feature_PA;
        case BaseFeature::PKMID2           : return feature_PI;
        case BaseFeature::PKMEAN2          : return feature_PS;
        case BaseFeature::LABEL_QV         : return feature_PQ;
        case BaseFeature::ALT_LABEL        : return feature_PT;
        case BaseFeature::ALT_LABEL_QV     : return feature_PV;
        case BaseFeature::PULSE_MERGE_QV   : return feature_PG;
        case BaseFeature::PULSE_CALL       : return feature_PC;
        case BaseFeature::PRE_PULSE_FRAMES : return feature_PD;
        case BaseFeature::PULSE_CALL_WIDTH : return feature_PX;
        case BaseFeature::START_FRAME      : return feature_SF;
        case BaseFeature::PULSE_EXCLUSION  : return feature_PE;
        default:
            throw std::runtime_error{ "unrecognized base feature" };
    }
}

static std::string FrameCodecName(const FrameCodec& codec)
{
    switch (codec) {
        case FrameCodec::RAW : return codec_RAW;
        case FrameCodec::V1  : return codec_V1;
        default:
            throw std::runtime_error{ "unrecognized frame codec" };
    }
}

static std::string BarcodeModeName(const BarcodeModeType& mode)
{
    switch (mode) {
        case BarcodeModeType::NONE       : return barcodemode_NONE;
        case BarcodeModeType::SYMMETRIC  : return barcodemode_SYM;
        case BarcodeModeType::ASYMMETRIC : return barcodemode_ASYM;
        case BarcodeModeType::TAILED     : return barcodemode_TAIL;
        default:
            throw std::runtime_error{ "unrecognized barcode mode type" };
    }
}

static std::string BarcodeQualityName(const BarcodeQualityType& type)
{
    switch (type) {
        case BarcodeQualityType::NONE        : return barcodequal_NONE;
        case BarcodeQualityType::SCORE       : return barcodequal_SCORE;
        case BarcodeQualityType::PROBABILITY : return barcodequal_PROB;
        default:
            throw std::runtime_error{ "unrecognized barcode quality type" };
    }
}

static std::string PlatformModelName(const PlatformModelType& type)
{
    switch (type) {
        case PlatformModelType::ASTRO  : return platformModelType_ASTRO;
        case PlatformModelType::RS     : return platformModelType_RS;
        case PlatformModelType::SEQUEL : return platformModelType_SEQUEL;
        default:
            throw std::runtime_error{ "unrecognized platform model type" };
    }
}

static const std::map<std::string, BaseFeature> nameToFeature
{
    { feature_DQ, BaseFeature::DELETION_QV },
    { feature_DT, BaseFeature::DELETION_TAG },
    { feature_IQ, BaseFeature::INSERTION_QV },
    { feature_MQ, BaseFeature::MERGE_QV },
    { feature_SQ, BaseFeature::SUBSTITUTION_QV },
    { feature_ST, BaseFeature::SUBSTITUTION_TAG },
    { feature_IP, BaseFeature::IPD },
    { feature_PW, BaseFeature::PULSE_WIDTH },
    { feature_PM, BaseFeature::PKMID },
    { feature_PA, BaseFeature::PKMEAN },
    { feature_PI, BaseFeature::PKMID2 },
    { feature_PS, BaseFeature::PKMEAN2 },
    { feature_PQ, BaseFeature::LABEL_QV },
    { feature_PT, BaseFeature::ALT_LABEL },
    { feature_PV, BaseFeature::ALT_LABEL_QV },
    { feature_PC, BaseFeature::PULSE_CALL },
    { feature_PG, BaseFeature::PULSE_MERGE_QV },
    { feature_PD, BaseFeature::PRE_PULSE_FRAMES },
    { feature_PX, BaseFeature::PULSE_CALL_WIDTH },
    { feature_SF, BaseFeature::START_FRAME },
    { feature_PE, BaseFeature::PULSE_EXCLUSION }
};

static const std::map<std::string, FrameCodec> nameToCodec
{
    { codec_RAW, FrameCodec::RAW },
    { codec_V1,  FrameCodec::V1 }
};

static const std::map<std::string, BarcodeModeType> nameToBarcodeMode
{
    { barcodemode_NONE, BarcodeModeType::NONE },
    { barcodemode_SYM,  BarcodeModeType::SYMMETRIC },
    { barcodemode_ASYM, BarcodeModeType::ASYMMETRIC },
    { barcodemode_TAIL, BarcodeModeType::TAILED }
};

static const std::map<std::string, BarcodeQualityType> nameToBarcodeQuality
{
    { barcodequal_NONE,  BarcodeQualityType::NONE },
    { barcodequal_SCORE, BarcodeQualityType::SCORE },
    { barcodequal_PROB,  BarcodeQualityType::PROBABILITY }
};

static const std::map<std::string, PlatformModelType> nameToPlatformModel
{
    { platformModelType_ASTRO,  PlatformModelType::ASTRO },
    { platformModelType_RS,     PlatformModelType::RS },
    { platformModelType_SEQUEL, PlatformModelType::SEQUEL }
};
// clang-format on

static inline bool IsLikelyBarcodeKey(const std::string& name) { return name.find("Barcode") == 0; }

static inline bool IsBaseFeature(const std::string& name)
{
    return nameToFeature.find(name) != nameToFeature.cend();
}

static inline BaseFeature BaseFeatureFromName(const std::string& name)
{
    return nameToFeature.at(name);
}

static inline FrameCodec FrameCodecFromName(const std::string& name)
{
    return nameToCodec.at(name);
}

static inline BarcodeModeType BarcodeModeFromName(const std::string& name)
{
    return nameToBarcodeMode.at(name);
}

static inline BarcodeQualityType BarcodeQualityFromName(const std::string& name)
{
    return nameToBarcodeQuality.at(name);
}

static inline PlatformModelType PlatformModelFromName(std::string name)
{
    return nameToPlatformModel.at(name);
}

}  // namespace internal

ReadGroupInfo::ReadGroupInfo() : readType_{"UNKNOWN"} {}

ReadGroupInfo::ReadGroupInfo(std::string id) : id_{std::move(id)}, readType_{"UNKNOWN"} {}

ReadGroupInfo::ReadGroupInfo(std::string movieName, std::string readType)
    : id_{MakeReadGroupId(movieName, readType)}
    , movieName_{std::move(movieName)}
    , readType_{std::move(readType)}
{
}

ReadGroupInfo::ReadGroupInfo(std::string movieName, std::string readType,
                             PlatformModelType platform)
    : id_{MakeReadGroupId(movieName, readType)}
    , movieName_{std::move(movieName)}
    , platformModel_{std::move(platform)}
    , readType_{std::move(readType)}
{
}

void ReadGroupInfo::DecodeSamDescription(const std::string& description)
{
    // split on semicolons
    // for each, split on equal
    //    determine name ->

    const auto tokens = internal::Split(description, ';');
    if (tokens.empty()) return;

    bool hasBarcodeFile = false;
    bool hasBarcodeHash = false;
    bool hasBarcodeCount = false;
    bool hasBarcodeMode = false;
    bool hasBarcodeQuality = false;

    // iterate over tokens
    for (const auto& token : tokens) {

        const auto foundEqual = token.find('=');
        if (foundEqual == std::string::npos) continue;

        const auto key = token.substr(0, foundEqual);
        auto value = token.substr(foundEqual + 1);

        // 'mandatory' items
        // clang-format off
        if      (key == internal::token_RT) readType_ = std::move(value);
        else if (key == internal::token_BK) bindingKit_ = std::move(value);
        else if (key == internal::token_BV) basecallerVersion_ = std::move(value);
        else if (key == internal::token_SK) sequencingKit_ = std::move(value);
        else if (key == internal::token_FR) frameRateHz_ = std::move(value);
        else if (key == internal::token_CT) control_ = (value == "TRUE");
        // clang-format on

        // base features
        else if (internal::IsBaseFeature(key))
            features_[internal::BaseFeatureFromName(key)] = std::move(value);

        // barcode data
        else if (internal::IsLikelyBarcodeKey(key)) {
            if (key == internal::token_BF) {
                barcodeFile_ = std::move(value);
                hasBarcodeFile = true;
            } else if (key == internal::token_BH) {
                barcodeHash_ = std::move(value);
                hasBarcodeHash = true;
            } else if (key == internal::token_BC) {
                barcodeCount_ = std::stoul(value);
                hasBarcodeCount = true;
            } else if (key == internal::token_BM) {
                barcodeMode_ = internal::BarcodeModeFromName(value);
                hasBarcodeMode = true;
            } else if (key == internal::token_BQ) {
                barcodeQuality_ = internal::BarcodeQualityFromName(value);
                hasBarcodeQuality = true;
            }
        }

        // frame codecs
        else {
            const auto keyParts = internal::Split(key, ':');
            if (keyParts.size() == 2) {
                const auto& subkey = keyParts.at(0);
                if (subkey == internal::feature_IP) {
                    ipdCodec_ = internal::FrameCodecFromName(keyParts.at(1));
                    features_[BaseFeature::IPD] = std::move(value);
                } else if (subkey == internal::feature_PW) {
                    pulseWidthCodec_ = internal::FrameCodecFromName(keyParts.at(1));
                    features_[BaseFeature::PULSE_WIDTH] = std::move(value);
                }
            }
        }
    }

    hasBarcodeData_ = (hasBarcodeFile && hasBarcodeHash && hasBarcodeCount && hasBarcodeMode &&
                       hasBarcodeQuality);
}

std::string ReadGroupInfo::EncodeSamDescription() const
{
    constexpr static const char SEP = ';';
    constexpr static const char COLON = ':';
    constexpr static const char EQ = '=';

    std::string result{internal::token_RT + EQ + readType_};

    std::string featureName;
    for (const auto& feature : features_) {

        featureName = internal::BaseFeatureName(feature.first);
        if (featureName.empty() || feature.second.empty())
            continue;
        else if (featureName == internal::feature_IP) {
            featureName.push_back(COLON);
            featureName.append(internal::FrameCodecName(ipdCodec_));
        } else if (featureName == internal::feature_PW) {
            featureName.push_back(COLON);
            featureName.append(internal::FrameCodecName(pulseWidthCodec_));
        }
        result.append(SEP + featureName + EQ + feature.second);
    }

    // clang-format off
    if (!bindingKit_.empty())        result.append(SEP + internal::token_BK + EQ + bindingKit_);
    if (!sequencingKit_.empty())     result.append(SEP + internal::token_SK + EQ + sequencingKit_);
    if (!basecallerVersion_.empty()) result.append(SEP + internal::token_BV + EQ + basecallerVersion_);
    if (!frameRateHz_.empty())       result.append(SEP + internal::token_FR + EQ + frameRateHz_);
    if (control_)                    result.append(SEP + internal::token_CT + EQ + (control_ ? "TRUE" : "FALSE"));
    // clang-format on

    if (hasBarcodeData_) {
        const std::string barcodeData{
            SEP + internal::token_BF + EQ + barcodeFile_ + SEP + internal::token_BH + EQ +
            barcodeHash_ + SEP + internal::token_BC + EQ + std::to_string(barcodeCount_) + SEP +
            internal::token_BM + EQ + internal::BarcodeModeName(barcodeMode_) + SEP +
            internal::token_BQ + EQ + internal::BarcodeQualityName(barcodeQuality_)};
        result.append(barcodeData);
    }

    return result;
}

ReadGroupInfo ReadGroupInfo::FromSam(const std::string& sam)
{
    // pop off '@RG\t', then split rest of line into tokens
    const auto tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty()) return {};

    ReadGroupInfo rg;
    std::map<std::string, std::string> custom;

    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set read group info
        // clang-format off
        if      (tokenTag == internal::sam_ID) rg.Id(std::move(tokenValue));
        else if (tokenTag == internal::sam_CN) rg.SequencingCenter(std::move(tokenValue));
        else if (tokenTag == internal::sam_DT) rg.Date(std::move(tokenValue));
        else if (tokenTag == internal::sam_FO) rg.FlowOrder(std::move(tokenValue));
        else if (tokenTag == internal::sam_KS) rg.KeySequence(std::move(tokenValue));
        else if (tokenTag == internal::sam_LB) rg.Library(std::move(tokenValue));
        else if (tokenTag == internal::sam_PG) rg.Programs(std::move(tokenValue));
        else if (tokenTag == internal::sam_PI) rg.PredictedInsertSize(std::move(tokenValue));
        else if (tokenTag == internal::sam_PU) rg.MovieName(std::move(tokenValue));
        else if (tokenTag == internal::sam_SM) rg.Sample(std::move(tokenValue));
        else if (tokenTag == internal::sam_DS) rg.DecodeSamDescription(std::move(tokenValue));
        else if (tokenTag == internal::sam_PM) rg.PlatformModel(internal::PlatformModelFromName(std::move(tokenValue)));
        // clang-format on

        // if not platform name (always "PACBIO" for us), store as a custom tag
        else if (tokenTag != internal::sam_PL)
            custom[tokenTag] = std::move(tokenValue);
    }
    rg.CustomTags(std::move(custom));

    return rg;
}

std::string ReadGroupInfo::IntToId(const int32_t id)
{
    std::ostringstream s;
    s << std::setfill('0') << std::setw(8) << std::hex << id;
    return s.str();
}

ReadGroupInfo& ReadGroupInfo::IpdCodec(FrameCodec codec, std::string tag)
{
    // store desired codec type
    ipdCodec_ = std::move(codec);

    // update base features map
    const std::string actualTag = (tag.empty() ? "ip" : std::move(tag));
    BaseFeatureTag(BaseFeature::IPD, actualTag);
    return *this;
}

ReadGroupInfo& ReadGroupInfo::PulseWidthCodec(FrameCodec codec, std::string tag)
{
    // store desired codec type
    pulseWidthCodec_ = std::move(codec);

    // update base features map
    const std::string actualTag = (tag.empty() ? "pw" : std::move(tag));
    BaseFeatureTag(BaseFeature::PULSE_WIDTH, actualTag);
    return *this;
}

std::string ReadGroupInfo::SequencingChemistryFromTriple(const std::string& bindingKit,
                                                         const std::string& sequencingKit,
                                                         const std::string& basecallerVersion)
{
    const auto verFields = internal::Split(basecallerVersion, '.');
    if (verFields.size() < 2)
        throw std::runtime_error{"basecaller version too short: " + basecallerVersion};
    const std::string version{verFields.at(0) + '.' + verFields.at(1)};

    // check updated table first, if it exists (empty if not), overriding the built-in lookup
    for (const auto& row : internal::GetChemistryTableFromEnv()) {
        if (bindingKit == row[0] && sequencingKit == row[1] && version == row[2]) return row[3];
    }

    for (const auto& row : internal::BuiltInChemistryTable) {
        if (bindingKit == row[0] && sequencingKit == row[1] && version == row[2]) return row[3];
    }

    // not found
    throw InvalidSequencingChemistryException{bindingKit, sequencingKit, basecallerVersion};
}

std::string ReadGroupInfo::ToSam() const
{
    std::ostringstream out;
    out << "@RG" << internal::MakeSamTag(internal::sam_ID, id_)
        << internal::MakeSamTag(internal::sam_PL, Platform());

    const auto description = EncodeSamDescription();
    if (!description.empty()) out << internal::MakeSamTag(internal::sam_DS, description);

    // clang-format off
    if (!sequencingCenter_.empty())    out << internal::MakeSamTag(internal::sam_CN, sequencingCenter_);
    if (!date_.empty())                out << internal::MakeSamTag(internal::sam_DT, date_);
    if (!flowOrder_.empty())           out << internal::MakeSamTag(internal::sam_FO, flowOrder_);
    if (!keySequence_.empty())         out << internal::MakeSamTag(internal::sam_KS, keySequence_);
    if (!library_.empty())             out << internal::MakeSamTag(internal::sam_LB, library_);
    if (!programs_.empty())            out << internal::MakeSamTag(internal::sam_PG, programs_);
    if (!predictedInsertSize_.empty()) out << internal::MakeSamTag(internal::sam_PI, predictedInsertSize_);
    if (!movieName_.empty())           out << internal::MakeSamTag(internal::sam_PU, movieName_);
    if (!sample_.empty())              out << internal::MakeSamTag(internal::sam_SM, sample_);
    // clang-format on

    out << internal::MakeSamTag(internal::sam_PM, internal::PlatformModelName(platformModel_));

    // append any custom tags
    for (const auto& attribute : custom_)
        out << internal::MakeSamTag(attribute.first, attribute.second);

    return out.str();
}

std::string MakeReadGroupId(const std::string& movieName, const std::string& readType)
{
    return MD5Hash(movieName + "//" + readType).substr(0, 8);
}

bool ReadGroupInfo::operator==(const ReadGroupInfo& other) const
{
    return id_ == other.id_ && sequencingCenter_ == other.sequencingCenter_ &&
           date_ == other.date_ && flowOrder_ == other.flowOrder_ &&
           keySequence_ == other.keySequence_ && library_ == other.library_ &&
           programs_ == other.programs_ && platformModel_ == other.platformModel_ &&
           predictedInsertSize_ == other.predictedInsertSize_ && movieName_ == other.movieName_ &&
           sample_ == other.sample_ && readType_ == other.readType_ &&
           bindingKit_ == other.bindingKit_ && sequencingKit_ == other.sequencingKit_ &&
           basecallerVersion_ == other.basecallerVersion_ && frameRateHz_ == other.frameRateHz_ &&
           control_ == other.control_ && ipdCodec_ == other.ipdCodec_ &&
           pulseWidthCodec_ == other.pulseWidthCodec_ && hasBarcodeData_ == other.hasBarcodeData_ &&
           barcodeFile_ == other.barcodeFile_ && barcodeHash_ == other.barcodeHash_ &&
           barcodeCount_ == other.barcodeCount_ && barcodeMode_ == other.barcodeMode_ &&
           barcodeQuality_ == other.barcodeQuality_ && features_.size() == other.features_.size() &&
           std::equal(features_.cbegin(), features_.cend(), other.features_.cbegin()) &&
           custom_.size() == other.custom_.size() &&
           std::equal(custom_.begin(), custom_.end(), other.custom_.cbegin());
}

}  // namespace BAM
}  // namespace PacBio

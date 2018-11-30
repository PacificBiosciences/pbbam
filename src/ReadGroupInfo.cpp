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
#include <tuple>
#include <unordered_map>

#include <boost/algorithm/cxx14/equal.hpp>
#include <boost/algorithm/string.hpp>

#include "ChemistryTable.h"
#include "EnumClassHash.h"
#include "pbbam/MD5.h"
#include "pbbam/SamTagCodec.h"
#include "pbbam/StringUtilities.h"

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

// clang-format off
static std::string BaseFeatureName(const BaseFeature& feature)
{
    static const std::unordered_map<BaseFeature, std::string, EnumClassHash> lookup{
        {BaseFeature::DELETION_QV,      feature_DQ},
        {BaseFeature::DELETION_TAG,     feature_DT},
        {BaseFeature::INSERTION_QV,     feature_IQ},
        {BaseFeature::MERGE_QV,         feature_MQ},
        {BaseFeature::SUBSTITUTION_QV,  feature_SQ},
        {BaseFeature::SUBSTITUTION_TAG, feature_ST},
        {BaseFeature::IPD,              feature_IP},
        {BaseFeature::PULSE_WIDTH,      feature_PW},
        {BaseFeature::PKMID,            feature_PM},
        {BaseFeature::PKMEAN,           feature_PA},
        {BaseFeature::PKMID2,           feature_PI},
        {BaseFeature::PKMEAN2,          feature_PS},
        {BaseFeature::LABEL_QV,         feature_PQ},
        {BaseFeature::ALT_LABEL,        feature_PT},
        {BaseFeature::ALT_LABEL_QV,     feature_PV},
        {BaseFeature::PULSE_MERGE_QV,   feature_PG},
        {BaseFeature::PULSE_CALL,       feature_PC},
        {BaseFeature::PRE_PULSE_FRAMES, feature_PD},
        {BaseFeature::PULSE_CALL_WIDTH, feature_PX},
        {BaseFeature::START_FRAME,      feature_SF},
        {BaseFeature::PULSE_EXCLUSION,  feature_PE}
    };

    const auto found = lookup.find(feature);
    if (found != lookup.cend())
        return found->second;
    throw std::runtime_error{ "unrecognized base feature" };
}

static std::string FrameCodecName(const FrameCodec& codec)
{
    static const std::unordered_map<FrameCodec, std::string, EnumClassHash> lookup{
        {FrameCodec::RAW, codec_RAW},
        {FrameCodec::V1,  codec_V1}
    };

    const auto found = lookup.find(codec);
    if (found != lookup.cend())
        return found->second;
    throw std::runtime_error{ "unrecognized frame codec" };
}

static std::string BarcodeModeName(const BarcodeModeType& mode)
{
    static const std::unordered_map<BarcodeModeType, std::string, EnumClassHash> lookup{
        {BarcodeModeType::NONE,       barcodemode_NONE},
        {BarcodeModeType::SYMMETRIC,  barcodemode_SYM},
        {BarcodeModeType::ASYMMETRIC, barcodemode_ASYM},
        {BarcodeModeType::TAILED,     barcodemode_TAIL}
    };

    const auto found = lookup.find(mode);
    if (found != lookup.cend())
        return found->second;
    throw std::runtime_error{ "unrecognized barcode mode type" };
}

static std::string BarcodeQualityName(const BarcodeQualityType& type)
{
    static const std::unordered_map<BarcodeQualityType, std::string, EnumClassHash> lookup{
        {BarcodeQualityType::NONE,        barcodequal_NONE},
        {BarcodeQualityType::SCORE,       barcodequal_SCORE},
        {BarcodeQualityType::PROBABILITY, barcodequal_PROB}
    };

    const auto found = lookup.find(type);
    if (found != lookup.cend())
        return found->second;
    throw std::runtime_error{ "unrecognized barcode quality type" };
}

static std::string PlatformModelName(const PlatformModelType& type)
{
    static const std::unordered_map<PlatformModelType, std::string, EnumClassHash> lookup{
        {PlatformModelType::ASTRO,    platformModelType_ASTRO},
        {PlatformModelType::RS,       platformModelType_RS},
        {PlatformModelType::SEQUEL,   platformModelType_SEQUEL},
        {PlatformModelType::SEQUELII, platformModelType_SEQUELII}
    };

    const auto found = lookup.find(type);
    if (found != lookup.cend())
        return found->second;
    throw std::runtime_error{ "unrecognized platform model type" };
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
    { platformModelType_ASTRO,    PlatformModelType::ASTRO },
    { platformModelType_RS,       PlatformModelType::RS },
    { platformModelType_SEQUEL,   PlatformModelType::SEQUEL },
    { platformModelType_SEQUELII, PlatformModelType::SEQUELII }
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

}  // anonymous

ReadGroupInfo::ReadGroupInfo(std::string baseId, std::pair<uint16_t, uint16_t> barcodes)

{
    std::ostringstream id;
    id << baseId << '/' << std::to_string(barcodes.first) << "--"
       << std::to_string(barcodes.second);
    id_ = id.str();
    baseId_ = std::move(baseId);
    barcodes_ = std::move(barcodes);
}

void ReadGroupInfo::DecodeBarcodeKey(const std::string& key, std::string value)
{
    if (key == token_BF)
        barcodeFile_ = std::move(value);
    else if (key == token_BH)
        barcodeHash_ = std::move(value);
    else if (key == token_BC)
        barcodeCount_ = std::stoul(value);
    else if (key == token_BM)
        barcodeMode_ = BarcodeModeFromName(value);
    else if (key == token_BQ)
        barcodeQuality_ = BarcodeQualityFromName(value);
}

void ReadGroupInfo::DecodeFrameCodecKey(const std::string& key, std::string value)
{
    const auto keyParts = Split(key, ':');
    if (keyParts.size() == 2) {
        const auto& subkey = keyParts.at(0);
        if (subkey == feature_IP) {
            ipdCodec_ = FrameCodecFromName(keyParts.at(1));
            features_[BaseFeature::IPD] = std::move(value);
        } else if (subkey == feature_PW) {
            pulseWidthCodec_ = FrameCodecFromName(keyParts.at(1));
            features_[BaseFeature::PULSE_WIDTH] = std::move(value);
        }
    }
}

void ReadGroupInfo::DecodeSamDescription(const std::string& description)
{
    const auto tokens = Split(description, ';');
    if (tokens.empty()) return;

    // iterate over tokens
    for (const auto& token : tokens) {

        const auto foundEqual = token.find('=');
        if (foundEqual == std::string::npos) continue;

        const auto key = token.substr(0, foundEqual);
        auto value = token.substr(foundEqual + 1);

        // 'mandatory' items
        // clang-format off
        if      (key == token_RT) readType_ = std::move(value);
        else if (key == token_BK) bindingKit_ = std::move(value);
        else if (key == token_BV) basecallerVersion_ = std::move(value);
        else if (key == token_SK) sequencingKit_ = std::move(value);
        else if (key == token_FR) frameRateHz_ = std::move(value);
        else if (key == token_CT) control_ = (value == "TRUE");
        // clang-format on

        // base features
        else if (IsBaseFeature(key))
            features_[BaseFeatureFromName(key)] = std::move(value);

        // barcode data
        else if (IsLikelyBarcodeKey(key))
            DecodeBarcodeKey(key, std::move(value));

        // frame codecs
        else
            DecodeFrameCodecKey(key, std::move(value));
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
        if (featureName.empty() || feature.second.empty())
            continue;
        else if (featureName == feature_IP) {
            featureName.push_back(COLON);
            featureName.append(FrameCodecName(ipdCodec_));
        } else if (featureName == feature_PW) {
            featureName.push_back(COLON);
            featureName.append(FrameCodecName(pulseWidthCodec_));
        }
        result.append(SEP + featureName + EQ + feature.second);
    }

    // clang-format off
    if (!bindingKit_.empty())        result.append(SEP + token_BK + EQ + bindingKit_);
    if (!sequencingKit_.empty())     result.append(SEP + token_SK + EQ + sequencingKit_);
    if (!basecallerVersion_.empty()) result.append(SEP + token_BV + EQ + basecallerVersion_);
    if (!frameRateHz_.empty())       result.append(SEP + token_FR + EQ + frameRateHz_);
    if (control_)                    result.append(SEP + token_CT + EQ + (control_ ? "TRUE" : "FALSE"));
    // clang-format on

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

ReadGroupInfo ReadGroupInfo::FromSam(const std::string& sam)
{
    // pop off '@RG\t', then split rest of line into tokens
    const auto tokens = Split(sam.substr(4), '\t');
    if (tokens.empty()) return {};

    ReadGroupInfo rg;
    std::map<std::string, std::string> custom;

    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set read group info
        // clang-format off
        if      (tokenTag == sam_ID) rg.Id(std::move(tokenValue));
        else if (tokenTag == sam_CN) rg.SequencingCenter(std::move(tokenValue));
        else if (tokenTag == sam_DT) rg.Date(std::move(tokenValue));
        else if (tokenTag == sam_FO) rg.FlowOrder(std::move(tokenValue));
        else if (tokenTag == sam_KS) rg.KeySequence(std::move(tokenValue));
        else if (tokenTag == sam_LB) rg.Library(std::move(tokenValue));
        else if (tokenTag == sam_PG) rg.Programs(std::move(tokenValue));
        else if (tokenTag == sam_PI) rg.PredictedInsertSize(std::move(tokenValue));
        else if (tokenTag == sam_PU) rg.MovieName(std::move(tokenValue));
        else if (tokenTag == sam_SM) rg.Sample(std::move(tokenValue));
        else if (tokenTag == sam_DS) rg.DecodeSamDescription(std::move(tokenValue));
        else if (tokenTag == sam_PM) rg.PlatformModel(PlatformModelFromName(std::move(tokenValue)));
        // clang-format on

        // if not platform name (always "PACBIO" for us), store as a custom tag
        else if (tokenTag != sam_PL)
            custom[tokenTag] = std::move(tokenValue);
    }
    rg.CustomTags(std::move(custom));

    return rg;
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
                "could not fetch barcodes from malformatted read group ID: " + id +
                " Must be in the form: {RGID_STRING}/{bcForward}--{bcReverse}"};
        }

        // catch here so we can give more informative message
        try {
            barcodes_ = std::pair<uint16_t, uint16_t>(static_cast<uint16_t>(std::stoul(tokens[0])),
                                                      static_cast<uint16_t>(std::stoul(tokens[2])));
        } catch (std::exception& e) {
            throw std::runtime_error{
                "could not fetch barcodes from malformatted read group ID: " + id_ +
                " Must be in the form: {RGID_STRING}/{bcForward}--{bcReverse}"};
        }
    }

    baseId_ = id.substr(0, slashAt);
    id_ = std::move(id);
    return *this;
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
    const auto verFields = Split(basecallerVersion, '.');
    if (verFields.size() < 2)
        throw std::runtime_error{"basecaller version too short: " + basecallerVersion};
    const std::string version{verFields.at(0) + '.' + verFields.at(1)};

    // check updated table first, if it exists (empty if not), overriding the built-in lookup
    for (const auto& row : GetChemistryTableFromEnv()) {
        if (bindingKit == row[0] && sequencingKit == row[1] && version == row[2]) return row[3];
    }

    for (const auto& row : BuiltInChemistryTable()) {
        if (bindingKit == row[0] && sequencingKit == row[1] && version == row[2]) return row[3];
    }

    // not found
    throw InvalidSequencingChemistryException{bindingKit, sequencingKit, basecallerVersion};
}

std::string ReadGroupInfo::ToSam() const
{
    std::ostringstream out;
    out << "@RG" << MakeSamTag(sam_ID, id_) << MakeSamTag(sam_PL, Platform());

    const auto description = EncodeSamDescription();
    if (!description.empty()) out << MakeSamTag(sam_DS, description);

    // clang-format off
    if (!sequencingCenter_.empty())    out << MakeSamTag(sam_CN, sequencingCenter_);
    if (!date_.empty())                out << MakeSamTag(sam_DT, date_);
    if (!flowOrder_.empty())           out << MakeSamTag(sam_FO, flowOrder_);
    if (!keySequence_.empty())         out << MakeSamTag(sam_KS, keySequence_);
    if (!library_.empty())             out << MakeSamTag(sam_LB, library_);
    if (!programs_.empty())            out << MakeSamTag(sam_PG, programs_);
    if (!predictedInsertSize_.empty()) out << MakeSamTag(sam_PI, predictedInsertSize_);
    if (!movieName_.empty())           out << MakeSamTag(sam_PU, movieName_);
    if (!sample_.empty())              out << MakeSamTag(sam_SM, sample_);
    if (barcodes_)
    {
        out << '\t' << sam_BC << ':'
            << barcodes_->first << "--" << barcodes_->second;
    }
    // clang-format on

    out << MakeSamTag(sam_PM, PlatformModelName(platformModel_));

    // append any custom tags
    for (const auto& attribute : custom_)
        out << MakeSamTag(attribute.first, attribute.second);

    return out.str();
}

std::string MakeReadGroupId(const std::string& movieName, const std::string& readType)
{
    return MD5Hash(movieName + "//" + readType).substr(0, 8);
}

bool ReadGroupInfo::operator==(const ReadGroupInfo& other) const
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

}  // namespace BAM
}  // namespace PacBio

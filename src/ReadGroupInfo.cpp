// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// File Description
/// \file ReadGroupInfo.cpp
/// \brief Implements the ReadGroupInfo class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ReadGroupInfo.h"
#include "pbbam/MD5.h"
#include "ChemistryTable.h"
#include "SequenceUtils.h"
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {
namespace internal {

static const std::string sam_ID = std::string{ "ID" };
static const std::string sam_CN = std::string{ "CN" };
static const std::string sam_DS = std::string{ "DS" };
static const std::string sam_DT = std::string{ "DT" };
static const std::string sam_FO = std::string{ "FO" };
static const std::string sam_KS = std::string{ "KS" };
static const std::string sam_LB = std::string{ "LB" };
static const std::string sam_PG = std::string{ "PG" };
static const std::string sam_PI = std::string{ "PI" };
static const std::string sam_PL = std::string{ "PL" };
static const std::string sam_PM = std::string{ "PM" };
static const std::string sam_PU = std::string{ "PU" };
static const std::string sam_SM = std::string{ "SM" };

static const std::string feature_DQ = std::string{ "DeletionQV" };
static const std::string feature_DT = std::string{ "DeletionTag" };
static const std::string feature_IQ = std::string{ "InsertionQV" };
static const std::string feature_MQ = std::string{ "MergeQV" };
static const std::string feature_SQ = std::string{ "SubstitutionQV" };
static const std::string feature_ST = std::string{ "SubstitutionTag" };
static const std::string feature_IP = std::string{ "Ipd" };
static const std::string feature_PW = std::string{ "PulseWidth" };
static const std::string feature_PM = std::string{ "PkMid" };
static const std::string feature_PA = std::string{ "PkMean" };
static const std::string feature_PI = std::string{ "PkMid2" };
static const std::string feature_PS = std::string{ "PkMean2" };
static const std::string feature_LT = std::string{ "Label" };
static const std::string feature_PQ = std::string{ "LabelQV" };
static const std::string feature_PT = std::string{ "AltLabel" };
static const std::string feature_PV = std::string{ "AltLabelQV" };
static const std::string feature_PG = std::string{ "PulseMergeQV" };
static const std::string feature_PC = std::string{ "PulseCall" };
static const std::string feature_PD = std::string{ "PrePulseFrames" };
static const std::string feature_PX = std::string{ "PulseCallWidth" };
static const std::string feature_SF = std::string{ "StartFrame" };

static const std::string token_RT = std::string{ "READTYPE" };
static const std::string token_BK = std::string{ "BINDINGKIT" };
static const std::string token_SK = std::string{ "SEQUENCINGKIT" };
static const std::string token_BV = std::string{ "BASECALLERVERSION" };
static const std::string token_FR = std::string{ "FRAMERATEHZ" };
static const std::string token_CT = std::string{ "CONTROL" };

static const std::string token_BF = std::string{ "BarcodeFile" };
static const std::string token_BH = std::string{ "BarcodeHash" };
static const std::string token_BC = std::string{ "BarcodeCount" };
static const std::string token_BM = std::string{ "BarcodeMode" };
static const std::string token_BQ = std::string{ "BarcodeQuality" };

static const std::string codec_RAW = std::string{ "Frames" };
static const std::string codec_V1  = std::string{ "CodecV1" };

static const std::string barcodemode_NONE = std::string{ "None" };
static const std::string barcodemode_SYM  = std::string{ "Symmetric" };
static const std::string barcodemode_ASYM = std::string{ "Asymmetric" };
static const std::string barcodemode_TAIL = std::string{ "Tailed" };

static const std::string barcodequal_NONE  = std::string{ "None" };
static const std::string barcodequal_SCORE = std::string{ "Score" };
static const std::string barcodequal_PROB  = std::string{ "Probability" };

static const std::string platformModelType_ASTRO  = std::string{ "ASTRO" };
static const std::string platformModelType_RS     = std::string{ "RS" };
static const std::string platformModelType_SEQUEL = std::string{ "SEQUEL" };

static std::string BaseFeatureName(const BaseFeature& feature)
{
    switch(feature) {
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
        case BarcodeQualityType::NONE  : return barcodequal_NONE;
        case BarcodeQualityType::SCORE : return barcodequal_SCORE;
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

static const auto nameToFeature = std::map<std::string, BaseFeature>
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
    { feature_SF, BaseFeature::START_FRAME }
};

static const auto nameToCodec = std::map<std::string, FrameCodec>
{
    { codec_RAW, FrameCodec::RAW },
    { codec_V1,  FrameCodec::V1 }
};

static const auto nameToBarcodeMode = std::map<std::string, BarcodeModeType>
{
    { barcodemode_NONE, BarcodeModeType::NONE },
    { barcodemode_SYM,  BarcodeModeType::SYMMETRIC },
    { barcodemode_ASYM, BarcodeModeType::ASYMMETRIC },
    { barcodemode_TAIL, BarcodeModeType::TAILED }
};

static const auto nameToBarcodeQuality = std::map<std::string, BarcodeQualityType>
{
    { barcodequal_NONE,  BarcodeQualityType::NONE },
    { barcodequal_SCORE, BarcodeQualityType::SCORE },
    { barcodequal_PROB,  BarcodeQualityType::PROBABILITY }
};

static const auto nameToPlatformModel = std::map<std::string, PlatformModelType>
{
    { platformModelType_ASTRO,  PlatformModelType::ASTRO },
    { platformModelType_RS,     PlatformModelType::RS },
    { platformModelType_SEQUEL, PlatformModelType::SEQUEL }
};

static inline bool IsLikelyBarcodeKey(const std::string& name)
{
    return name.find("Barcode") == 0;
}

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

static inline PlatformModelType PlatformModelFromName(const std::string& name)
{
    return nameToPlatformModel.at(name);
}

} // namespace internal

ReadGroupInfo::ReadGroupInfo(void)
    : platformModel_(PlatformModelType::SEQUEL)
    , readType_("UNKNOWN")
    , ipdCodec_(FrameCodec::V1)
    , pulseWidthCodec_(FrameCodec::V1)
{ }

ReadGroupInfo::ReadGroupInfo(const std::string& id)
    : id_(id)
    , platformModel_(PlatformModelType::SEQUEL)
    , readType_("UNKNOWN")
    , ipdCodec_(FrameCodec::V1)
    , pulseWidthCodec_(FrameCodec::V1)
{ }

ReadGroupInfo::ReadGroupInfo(const std::string& movieName,
                             const std::string& readType)
    : id_(MakeReadGroupId(movieName, readType))
    , movieName_(movieName)
    , platformModel_(PlatformModelType::SEQUEL)
    , readType_(readType)
    , ipdCodec_(FrameCodec::V1)
    , pulseWidthCodec_(FrameCodec::V1)
{ }

ReadGroupInfo::ReadGroupInfo(const std::string& movieName,
                             const std::string& readType,
                             const PlatformModelType platform)
    : id_(MakeReadGroupId(movieName, readType))
    , movieName_(movieName)
    , platformModel_(platform)
    , readType_(readType)
    , ipdCodec_(FrameCodec::V1)
    , pulseWidthCodec_(FrameCodec::V1)
{ }

void ReadGroupInfo::DecodeSamDescription(const std::string& description)
{
    // split on semicolons
    // for each, split on equal
    //    determine name ->

    auto tokens = internal::Split(description, ';');
    if (tokens.empty())
        return;

    bool hasBarcodeFile = false;
    bool hasBarcodeHash = false;
    bool hasBarcodeCount = false;
    bool hasBarcodeMode = false;
    bool hasBarcodeQuality = false;

    // iterate over tokens
    for (auto&& token : tokens) {

        const auto foundEqual = token.find('=');
        if (foundEqual == std::string::npos)
            continue;

        const auto key = token.substr(0,foundEqual);
        const auto value = token.substr(foundEqual+1);

        // 'mandatory' items
        if      (key == internal::token_RT) readType_ = value;
        else if (key == internal::token_BK) bindingKit_ = value;
        else if (key == internal::token_BV) basecallerVersion_ = value;
        else if (key == internal::token_SK) sequencingKit_ = value;
        else if (key == internal::token_FR) frameRateHz_ = value;
        else if (key == internal::token_CT) control_ = (value == "TRUE");

        // base features
        else if (internal::IsBaseFeature(key))
            features_[internal::BaseFeatureFromName(key)] = value;

        // barcode data
        else if (internal::IsLikelyBarcodeKey(key)) {
            if (key == internal::token_BF) {
                barcodeFile_ = value;
                hasBarcodeFile = true;
            }
            else if (key == internal::token_BH) {
                barcodeHash_ = value;
                hasBarcodeHash = true;
            }
            else if (key == internal::token_BC) {
                barcodeCount_ = static_cast<size_t>(std::stoul(value));
                hasBarcodeCount = true;
            }
            else if (key == internal::token_BM) {
                barcodeMode_ = internal::BarcodeModeFromName(value);
                hasBarcodeMode = true;
            }
            else if (key == internal::token_BQ) {
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
                    features_[BaseFeature::IPD] = value;
                } 
                else if (subkey == internal::feature_PW) {
                    pulseWidthCodec_ = internal::FrameCodecFromName(keyParts.at(1));
                    features_[BaseFeature::PULSE_WIDTH] = value;
                }
            }
        }
    }

    hasBarcodeData_ = (hasBarcodeFile  &&
                       hasBarcodeHash  &&
                       hasBarcodeCount &&
                       hasBarcodeMode  &&
                       hasBarcodeQuality);
}

std::string ReadGroupInfo::EncodeSamDescription(void) const
{
    auto result = std::string{ };
    result.reserve(256);
    result.append(std::string(internal::token_RT+"=" + readType_));

    static const auto SEP   = std::string{";"};
    static const auto COLON = std::string{":"};
    static const auto EQ    = std::string{"="};

    auto featureName = std::string{ };
    const auto featureEnd = features_.cend();
    auto featureIter = features_.cbegin();
    for ( ; featureIter != featureEnd; ++featureIter ) {
        featureName = internal::BaseFeatureName(featureIter->first);
        if (featureName.empty() || featureIter->second.empty())
            continue;
        else if (featureName == internal::feature_IP) {
            featureName.append(COLON);
            featureName.append(internal::FrameCodecName(ipdCodec_));
        }
        else if (featureName == internal::feature_PW) {
            featureName.append(COLON);
            featureName.append(internal::FrameCodecName(pulseWidthCodec_));
        }
        result.append(std::string(SEP + featureName + EQ + featureIter->second));
    }

    if (!bindingKit_.empty())        result.append(SEP + internal::token_BK +EQ + bindingKit_);
    if (!sequencingKit_.empty())     result.append(SEP + internal::token_SK +EQ + sequencingKit_);
    if (!basecallerVersion_.empty()) result.append(SEP + internal::token_BV +EQ + basecallerVersion_);
    if (!frameRateHz_.empty())       result.append(SEP + internal::token_FR +EQ + frameRateHz_);
    if (control_)                    result.append(SEP + internal::token_CT +EQ + (control_ ? "TRUE"
                                                                                            : "FALSE"));

    if (hasBarcodeData_) {
        const auto barcodeData =
            std::string {
                SEP + internal::token_BF + EQ + barcodeFile_ +
                SEP + internal::token_BH + EQ + barcodeHash_ +
                SEP + internal::token_BC + EQ + std::to_string(barcodeCount_) +
                SEP + internal::token_BM + EQ + internal::BarcodeModeName(barcodeMode_) +
                SEP + internal::token_BQ + EQ + internal::BarcodeQualityName(barcodeQuality_)
            };
        result.reserve(result.size() + barcodeData.size());
        result.append(barcodeData);
    }

    return result;
}

ReadGroupInfo ReadGroupInfo::FromSam(const std::string& sam)
{
    // pop off '@RG\t', then split rest of line into tokens
    const auto tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty())
        return ReadGroupInfo{ };

    auto rg = ReadGroupInfo{ };
    auto custom = std::map<std::string, std::string>{ };

    for (auto&& token : tokens) {
        const auto tokenTag   = token.substr(0,2);
        const auto tokenValue = token.substr(3);

        // set read group info
        if      (tokenTag == internal::sam_ID) rg.Id(tokenValue);
        else if (tokenTag == internal::sam_CN) rg.SequencingCenter(tokenValue);
        else if (tokenTag == internal::sam_DT) rg.Date(tokenValue);
        else if (tokenTag == internal::sam_FO) rg.FlowOrder(tokenValue);
        else if (tokenTag == internal::sam_KS) rg.KeySequence(tokenValue);
        else if (tokenTag == internal::sam_LB) rg.Library(tokenValue);
        else if (tokenTag == internal::sam_PG) rg.Programs(tokenValue);
        else if (tokenTag == internal::sam_PI) rg.PredictedInsertSize(tokenValue);
        else if (tokenTag == internal::sam_PU) rg.MovieName(tokenValue);
        else if (tokenTag == internal::sam_SM) rg.Sample(tokenValue);
        else if (tokenTag == internal::sam_DS) rg.DecodeSamDescription(tokenValue);
        else if (tokenTag == internal::sam_PM) rg.PlatformModel(internal::PlatformModelFromName(tokenValue));

        // if not platform name (always "PACBIO" for us), store as a custom tag
        else if (tokenTag != internal::sam_PL)
            custom[tokenTag] = tokenValue;
    }
    rg.CustomTags(custom);

    return rg;
}

std::string ReadGroupInfo::IntToId(const int32_t id)
{
    std::stringstream s;
    s << std::setfill('0') << std::setw(8) << std::hex << id;
    return s.str();
}

ReadGroupInfo& ReadGroupInfo::IpdCodec(const FrameCodec& codec,
                                       const std::string& tag)
{
    // store desired codec type
    ipdCodec_ = codec;

    // update base features map
    auto actualTag = tag;
    if (actualTag.empty())
        actualTag = "ip";
    BaseFeatureTag(BaseFeature::IPD, actualTag);
    return *this;
}

ReadGroupInfo& ReadGroupInfo::PulseWidthCodec(const FrameCodec& codec,
                                              const std::string& tag)
{
    // store desired codec type
    pulseWidthCodec_ = codec;

    // update base features map
    auto actualTag = tag;
    if (actualTag.empty())
        actualTag = "pw";
    BaseFeatureTag(BaseFeature::PULSE_WIDTH, actualTag);
    return *this;
}

std::string ReadGroupInfo::SequencingChemistryFromTriple(const std::string& bindingKit,
                                                         const std::string& sequencingKit,
                                                         const std::string& basecallerVersion)
{
    const auto verFields = internal::Split(basecallerVersion, '.');
    if (verFields.size() < 2)
        throw std::runtime_error("basecaller version too short: " + basecallerVersion);
    const std::string ver = verFields.at(0) + "." + verFields.at(1);

    // check updated table first, if it exists (empty if not), overriding the built-in lookup
    for (const auto& row : internal::GetChemistryTableFromEnv()) {
        if (bindingKit == row[0] && sequencingKit == row[1] && ver == row[2])
            return row[3];
    }

    for (const auto& row : internal::BuiltInChemistryTable) {
        if (bindingKit == row[0] && sequencingKit == row[1] && ver == row[2])
            return row[3];
    }

    // not found
    throw InvalidSequencingChemistryException(bindingKit,
                                              sequencingKit,
                                              basecallerVersion);
}

std::string ReadGroupInfo::ToSam(void) const
{
    std::stringstream out;
    out << "@RG"
        << internal::MakeSamTag(internal::sam_ID, id_)
        << internal::MakeSamTag(internal::sam_PL, Platform());

    auto description = EncodeSamDescription();
    if (!description.empty())
        out << internal::MakeSamTag(internal::sam_DS, description);

    if (!sequencingCenter_.empty())    out << internal::MakeSamTag(internal::sam_CN, sequencingCenter_);
    if (!date_.empty())                out << internal::MakeSamTag(internal::sam_DT, date_);
    if (!flowOrder_.empty())           out << internal::MakeSamTag(internal::sam_FO, flowOrder_);
    if (!keySequence_.empty())         out << internal::MakeSamTag(internal::sam_KS, keySequence_);
    if (!library_.empty())             out << internal::MakeSamTag(internal::sam_LB, library_);
    if (!programs_.empty())            out << internal::MakeSamTag(internal::sam_PG, programs_);
    if (!predictedInsertSize_.empty()) out << internal::MakeSamTag(internal::sam_PI, predictedInsertSize_);
    if (!movieName_.empty())           out << internal::MakeSamTag(internal::sam_PU, movieName_);
    if (!sample_.empty())              out << internal::MakeSamTag(internal::sam_SM, sample_);

    out << internal::MakeSamTag(internal::sam_PM, internal::PlatformModelName(platformModel_));

    // append any custom tags
    auto customIter = custom_.cbegin();
    auto customEnd  = custom_.cend();
    for ( ; customIter != customEnd; ++customIter )
        out << internal::MakeSamTag(customIter->first, customIter->second);

    return out.str();
}

std::string MakeReadGroupId(const std::string& movieName,
                            const std::string& readType)
{
    return MD5Hash(movieName + "//" + readType).substr(0,8);
}

bool ReadGroupInfo::operator==(const ReadGroupInfo& other) const
{
    return id_ == other.id_ 
            && sequencingCenter_ == other.sequencingCenter_        
            && date_ == other.date_                    
            && flowOrder_ == other.flowOrder_               
            && keySequence_ == other.keySequence_             
            && library_ == other.library_                 
            && programs_ == other.programs_                
            && platformModel_ == other.platformModel_                
            && predictedInsertSize_ == other.predictedInsertSize_     
            && movieName_ == other.movieName_               
            && sample_ == other.sample_                  
            && readType_ == other.readType_ 
            && bindingKit_ == other.bindingKit_ 
            && sequencingKit_ == other.sequencingKit_ 
            && basecallerVersion_ == other.basecallerVersion_ 
            && frameRateHz_ == other.frameRateHz_ 
            && control_ == other.control_ 
            && ipdCodec_ == other.ipdCodec_
            && pulseWidthCodec_ == other.pulseWidthCodec_
            && hasBarcodeData_ == other.hasBarcodeData_
            && barcodeFile_ == other.barcodeFile_
            && barcodeHash_ == other.barcodeHash_
            && barcodeCount_ == other.barcodeCount_
            && barcodeMode_ == other.barcodeMode_
            && barcodeQuality_ == other.barcodeQuality_
            && features_.size() == other.features_.size()
            && std::equal(features_.cbegin(),
                          features_.cend(),
                          other.features_.cbegin())
            && custom_.size() == other.custom_.size()
            && std::equal(custom_.begin(),
                          custom_.end(),
                          other.custom_.cbegin());
}

} // namespace BAM
} // namespace PacBio

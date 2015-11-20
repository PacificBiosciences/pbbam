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

#include "pbbam/ReadGroupInfo.h"
#include "ChemistryTable.h"
#include "SequenceUtils.h"
#include <cram/md5.h>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <cstdio>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static const string sam_ID = string{ "ID" };
static const string sam_CN = string{ "CN" };
static const string sam_DS = string{ "DS" };
static const string sam_DT = string{ "DT" };
static const string sam_FO = string{ "FO" };
static const string sam_KS = string{ "KS" };
static const string sam_LB = string{ "LB" };
static const string sam_PG = string{ "PG" };
static const string sam_PI = string{ "PI" };
static const string sam_PL = string{ "PL" };
static const string sam_PU = string{ "PU" };
static const string sam_SM = string{ "SM" };

static const string feature_DQ = string{ "DeletionQV" };
static const string feature_DT = string{ "DeletionTag" };
static const string feature_IQ = string{ "InsertionQV" };
static const string feature_MQ = string{ "MergeQV" };
static const string feature_SQ = string{ "SubstitutionQV" };
static const string feature_ST = string{ "SubstitutionTag" };
static const string feature_IP = string{ "Ipd" };
static const string feature_PW = string{ "PulseWidth" };
static const string feature_PM = string{ "PkMid" };
static const string feature_PA = string{ "PkMean" };
static const string feature_LT = string{ "Label" };
static const string feature_PQ = string{ "LabelQV" };
static const string feature_PT = string{ "AltLabel" };
static const string feature_PV = string{ "AltLabelQV" };
static const string feature_PG = string{ "PulseMergeQV" };
static const string feature_PC = string{ "PulseCall" };
static const string feature_PD = string{ "PrePulseFrames" };
static const string feature_PX = string{ "PulseCallWidth" };
static const string feature_SF = string{ "StartFrame" };

static const string token_RT = string{ "READTYPE" };
static const string token_BK = string{ "BINDINGKIT" };
static const string token_SK = string{ "SEQUENCINGKIT" };
static const string token_BV = string{ "BASECALLERVERSION" };
static const string token_FR = string{ "FRAMERATEHZ" };
static const string token_CT = string{ "CONTROL" };

static const string token_BF = string{ "BarcodeFile" };
static const string token_BH = string{ "BarcodeHash" };
static const string token_BC = string{ "BarcodeCount" };
static const string token_BM = string{ "BarcodeMode" };
static const string token_BQ = string{ "BarcodeQuality" };

static const string codec_RAW = string{ "Frames" };
static const string codec_V1  = string{ "CodecV1" };

static const string barcodemode_NONE = string{ "None" };
static const string barcodemode_SYM  = string{ "Symmetric" };
static const string barcodemode_ASYM = string{ "Asymmetric" };

static const string barcodequal_NONE  = string{ "None" };
static const string barcodequal_SCORE = string{ "Score" };
static const string barcodequal_PROB  = string{ "Probability" };

static
string BaseFeatureName(const BaseFeature& feature)
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
        case BaseFeature::LABEL_QV         : return feature_PQ;
        case BaseFeature::ALT_LABEL        : return feature_PT;
        case BaseFeature::ALT_LABEL_QV     : return feature_PV;
        case BaseFeature::PULSE_MERGE_QV   : return feature_PG;
        case BaseFeature::PULSE_CALL       : return feature_PC;
        case BaseFeature::PRE_PULSE_FRAMES : return feature_PD;
        case BaseFeature::PULSE_CALL_WIDTH : return feature_PX;
    case BaseFeature::START_FRAME          : return feature_SF;
        default:
            throw std::runtime_error{ "unrecognized base feature" };
    }
    return string{ }; // unreachable
}

static
string FrameCodecName(const FrameCodec& codec)
{
    switch (codec) {
        case FrameCodec::RAW : return codec_RAW;
        case FrameCodec::V1  : return codec_V1;
        default:
            throw std::runtime_error{ "unrecognized frame codec" };
    }
    return string{ }; // unreachable
}

static
string BarcodeModeName(const BarcodeModeType& mode)
{
    switch (mode) {
        case BarcodeModeType::NONE       : return barcodemode_NONE;
        case BarcodeModeType::SYMMETRIC  : return barcodemode_SYM;
        case BarcodeModeType::ASYMMETRIC : return barcodemode_ASYM;
        default:
            throw std::runtime_error{ "unrecognized barcode mode" };
    }
    return string{ }; // unreachable
}

static
string BarcodeQualityName(const BarcodeQualityType& type)
{
    switch (type) {
        case BarcodeQualityType::NONE  : return barcodequal_NONE;
        case BarcodeQualityType::SCORE : return barcodequal_SCORE;
        case BarcodeQualityType::PROBABILITY : return barcodequal_PROB;
        default:
            throw std::runtime_error{ "unrecognized barcode quality type" };
    }
    return string{ }; // unreachable
}

static map<string, BaseFeature>        nameToFeature;
static map<string, FrameCodec>         nameToCodec;
static map<string, BarcodeModeType>    nameToBarcodeMode;
static map<string, BarcodeQualityType> nameToBarcodeQuality;

static inline
void InitNameToFeature(void)
{
    if (nameToFeature.empty()) {
        nameToFeature[feature_DQ] = BaseFeature::DELETION_QV;
        nameToFeature[feature_DT] = BaseFeature::DELETION_TAG;
        nameToFeature[feature_IQ] = BaseFeature::INSERTION_QV;
        nameToFeature[feature_MQ] = BaseFeature::MERGE_QV;
        nameToFeature[feature_SQ] = BaseFeature::SUBSTITUTION_QV;
        nameToFeature[feature_ST] = BaseFeature::SUBSTITUTION_TAG;
        nameToFeature[feature_IP] = BaseFeature::IPD;
        nameToFeature[feature_PW] = BaseFeature::PULSE_WIDTH;
        nameToFeature[feature_PM] = BaseFeature::PKMID;
        nameToFeature[feature_PA] = BaseFeature::PKMEAN;
        nameToFeature[feature_PQ] = BaseFeature::LABEL_QV;
        nameToFeature[feature_PT] = BaseFeature::ALT_LABEL;
        nameToFeature[feature_PV] = BaseFeature::ALT_LABEL_QV;
        nameToFeature[feature_PC] = BaseFeature::PULSE_CALL;
        nameToFeature[feature_PG] = BaseFeature::PULSE_MERGE_QV;
        nameToFeature[feature_PD] = BaseFeature::PRE_PULSE_FRAMES;
        nameToFeature[feature_PX] = BaseFeature::PULSE_CALL_WIDTH;
        nameToFeature[feature_SF] = BaseFeature::START_FRAME;
    }
}

static inline
void InitNameToCodec(void)
{
    if (nameToCodec.empty()) {
        nameToCodec[codec_RAW] = FrameCodec::RAW;
        nameToCodec[codec_V1]  = FrameCodec::V1;
    }
}

static inline
void InitNameToBarcodeMode(void)
{
    if (nameToBarcodeMode.empty()) {
        nameToBarcodeMode[barcodemode_NONE] = BarcodeModeType::NONE;
        nameToBarcodeMode[barcodemode_SYM]  = BarcodeModeType::SYMMETRIC;
        nameToBarcodeMode[barcodemode_ASYM] = BarcodeModeType::ASYMMETRIC;
    }
}

static inline
void InitNameToBarcodeQuality(void)
{
    if (nameToBarcodeQuality.empty()) {
        nameToBarcodeQuality[barcodequal_NONE]  = BarcodeQualityType::NONE;
        nameToBarcodeQuality[barcodequal_SCORE] = BarcodeQualityType::SCORE;
        nameToBarcodeQuality[barcodequal_PROB]  = BarcodeQualityType::PROBABILITY;
    }
}

static inline
bool IsLikelyBarcodeKey(const string& name)
{ return name.find("Barcode") == 0; }

static inline
bool IsBaseFeature(const string& name)
{
    InitNameToFeature();
    return nameToFeature.find(name) != nameToFeature.cend();
}

static inline
BaseFeature BaseFeatureFromName(const string& name)
{
    InitNameToFeature();
    return nameToFeature.at(name);
}

static inline
FrameCodec FrameCodecFromName(const string& name)
{
    InitNameToCodec();
    return nameToCodec.at(name);
}

static inline
BarcodeModeType BarcodeModeFromName(const string& name)
{
    InitNameToBarcodeMode();
    return nameToBarcodeMode.at(name);
}

static inline
BarcodeQualityType BarcodeQualityFromName(const string& name)
{
    InitNameToBarcodeQuality();
    return nameToBarcodeQuality.at(name);
}

} // namespace internal

ReadGroupInfo::ReadGroupInfo(void)
    : readType_("UNKNOWN")
    , ipdCodec_(FrameCodec::V1)
    , pulseWidthCodec_(FrameCodec::V1)
{ }

ReadGroupInfo::ReadGroupInfo(const std::string& id)
    : id_(id)
    , readType_("UNKNOWN")
    , ipdCodec_(FrameCodec::V1)
    , pulseWidthCodec_(FrameCodec::V1)
{ }

ReadGroupInfo::ReadGroupInfo(const std::string& movieName,
                             const std::string& readType)
    : id_(MakeReadGroupId(movieName, readType))
    , movieName_(movieName)
    , readType_(readType)
{ }

ReadGroupInfo::ReadGroupInfo(const ReadGroupInfo& other)
    : id_(other.id_)
    , sequencingCenter_(other.sequencingCenter_)
    , date_(other.date_)
    , flowOrder_(other.flowOrder_)
    , keySequence_(other.keySequence_)
    , library_(other.library_)
    , programs_(other.programs_)
    , predictedInsertSize_(other.predictedInsertSize_)
    , movieName_(other.movieName_)
    , sample_(other.sample_)
    , readType_(other.readType_)
    , bindingKit_(other.bindingKit_)
    , sequencingKit_(other.sequencingKit_)
    , basecallerVersion_(other.basecallerVersion_)
    , frameRateHz_(other.frameRateHz_)
    , control_(other.control_)
    , ipdCodec_(other.ipdCodec_)
    , pulseWidthCodec_(other.pulseWidthCodec_)
    , hasBarcodeData_(other.hasBarcodeData_)
    , barcodeFile_(other.barcodeFile_)
    , barcodeHash_(other.barcodeHash_)
    , barcodeCount_(other.barcodeCount_)
    , barcodeMode_(other.barcodeMode_)
    , barcodeQuality_(other.barcodeQuality_)
    , features_(other.features_)
{  }

ReadGroupInfo::ReadGroupInfo(ReadGroupInfo&& other)
    : id_(std::move(other.id_))
    , sequencingCenter_(std::move(other.sequencingCenter_))
    , date_(std::move(other.date_))
    , flowOrder_(std::move(other.flowOrder_))
    , keySequence_(std::move(other.keySequence_))
    , library_(std::move(other.library_))
    , programs_(std::move(other.programs_))
    , predictedInsertSize_(std::move(other.predictedInsertSize_))
    , movieName_(std::move(other.movieName_))
    , sample_(std::move(other.sample_))
    , readType_(std::move(other.readType_))
    , bindingKit_(std::move(other.bindingKit_))
    , sequencingKit_(std::move(other.sequencingKit_))
    , basecallerVersion_(std::move(other.basecallerVersion_))
    , frameRateHz_(std::move(other.frameRateHz_))
    , control_(std::move(other.control_))
    , ipdCodec_(std::move(other.ipdCodec_))
    , pulseWidthCodec_(std::move(other.pulseWidthCodec_))
    , hasBarcodeData_(std::move(other.hasBarcodeData_))
    , barcodeFile_(std::move(other.barcodeFile_))
    , barcodeHash_(std::move(other.barcodeHash_))
    , barcodeCount_(std::move(other.barcodeCount_))
    , barcodeMode_(std::move(other.barcodeMode_))
    , barcodeQuality_(std::move(other.barcodeQuality_))
    , features_(std::move(other.features_))
{ }

ReadGroupInfo::~ReadGroupInfo(void) { }

ReadGroupInfo& ReadGroupInfo::operator=(const ReadGroupInfo& other)
{
    id_ = other.id_;
    sequencingCenter_ = other.sequencingCenter_;
    date_ = other.date_;
    flowOrder_ = other.flowOrder_;
    keySequence_ = other.keySequence_;
    library_ = other.library_;
    programs_ = other.programs_;
    predictedInsertSize_ = other.predictedInsertSize_;
    movieName_ = other.movieName_;
    sample_ = other.sample_;
    readType_ = other.readType_;
    bindingKit_ = other.bindingKit_;
    sequencingKit_ = other.sequencingKit_;
    basecallerVersion_ = other.basecallerVersion_;
    frameRateHz_ = other.frameRateHz_;
    control_ = other.control_;
    ipdCodec_ = other.ipdCodec_;
    pulseWidthCodec_ = other.pulseWidthCodec_;
    hasBarcodeData_ = other.hasBarcodeData_;
    barcodeFile_  = other.barcodeFile_;
    barcodeHash_ = other.barcodeHash_;
    barcodeCount_ = other.barcodeCount_;
    barcodeMode_ = other.barcodeMode_;
    barcodeQuality_ = other.barcodeQuality_;
    features_ = other.features_;
    return *this;
}

ReadGroupInfo& ReadGroupInfo::operator=(ReadGroupInfo&& other)
{
    id_ = std::move(other.id_);
    sequencingCenter_ = std::move(other.sequencingCenter_);
    date_ = std::move(other.date_);
    flowOrder_ = std::move(other.flowOrder_);
    keySequence_ = std::move(other.keySequence_);
    library_ = std::move(other.library_);
    programs_ = std::move(other.programs_);
    predictedInsertSize_ = std::move(other.predictedInsertSize_);
    movieName_ = std::move(other.movieName_);
    sample_ = std::move(other.sample_);
    readType_ = std::move(other.readType_);
    bindingKit_ = std::move(other.bindingKit_);
    sequencingKit_ = std::move(other.sequencingKit_);
    basecallerVersion_ = std::move(other.basecallerVersion_);
    frameRateHz_ = std::move(other.frameRateHz_);
    control_ = std::move(other.control_);
    ipdCodec_ = std::move(other.ipdCodec_);
    pulseWidthCodec_ = std::move(other.pulseWidthCodec_);
    hasBarcodeData_ = std::move(other.hasBarcodeData_);
    barcodeFile_  = std::move(other.barcodeFile_);
    barcodeHash_ = std::move(other.barcodeHash_);
    barcodeCount_ = std::move(other.barcodeCount_);
    barcodeMode_ = std::move(other.barcodeMode_);
    barcodeQuality_ = std::move(other.barcodeQuality_);
    features_ = std::move(other.features_);
    return *this;
}

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
        if (foundEqual == string::npos)
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
    auto result = string{ };
    result.reserve(256);
    result.append(std::string(internal::token_RT+"=" + readType_));

    static const auto SEP   = string{";"};
    static const auto COLON = string{":"};
    static const auto EQ    = string{"="};

    auto featureName = string{ };
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
        result.append(string(SEP + featureName + EQ + featureIter->second));
    }

    if (!bindingKit_.empty())        result.append(SEP + internal::token_BK +EQ + bindingKit_);
    if (!sequencingKit_.empty())     result.append(SEP + internal::token_SK +EQ + sequencingKit_);
    if (!basecallerVersion_.empty()) result.append(SEP + internal::token_BV +EQ + basecallerVersion_);
    if (!frameRateHz_.empty())       result.append(SEP + internal::token_FR +EQ + frameRateHz_);
    if (control_)                    result.append(SEP + internal::token_CT +EQ + (control_ ? "TRUE"
                                                                                            : "FALSE"));

    if (hasBarcodeData_) {
        const auto barcodeData =
            string {
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

ReadGroupInfo ReadGroupInfo::FromSam(const string& sam)
{
    // pop off '@RG\t', then split rest of line into tokens
    const auto tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty())
        return ReadGroupInfo{ };

    auto rg = ReadGroupInfo{ };
    auto custom = map<string, string>{ };

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

        // otherwise, "custom" tag
        else
            custom[tokenTag] = tokenValue;
    }
    rg.CustomTags(custom);

    return rg;
}

string ReadGroupInfo::IntToId(const int32_t id)
{
    stringstream s;
    s << std::setfill('0') << std::setw(8) << std::hex << id;
    return s.str();
}

ReadGroupInfo& ReadGroupInfo::IpdCodec(const FrameCodec& codec,
                                       const string& tag)
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
                                              const string& tag)
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

string ReadGroupInfo::SequencingChemistryFromTriple(const string& bindingKit,
                                                    const string& sequencingKit,
                                                    const string& basecallerVersion)
{
    const string ver{ basecallerVersion.substr(0, 3) };
    for (const auto& row : internal::ChemistryTable) {
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
    stringstream out;
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
    MD5_CTX md5;
    unsigned char digest[16];
    char hexdigest[9];

    MD5_Init(&md5);
    MD5_Update(&md5, reinterpret_cast<void*>(const_cast<char*>(movieName.c_str())), movieName.size());
    MD5_Update(&md5, reinterpret_cast<void*>(const_cast<char*>("//")), 2);
    MD5_Update(&md5, reinterpret_cast<void*>(const_cast<char*>(readType.c_str())), readType.size());
    MD5_Final(digest, &md5);

    for (int i = 0; i < 4; ++i)
        sprintf(&hexdigest[2*i], "%02x", digest[i]);

    return std::string{hexdigest, 8};
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

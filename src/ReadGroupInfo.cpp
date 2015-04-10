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

// Author: Derek Barnett

#include "pbbam/ReadGroupInfo.h"
#include "SequenceUtils.h"
#include <set>
#include <sstream>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace PacBio {
namespace BAM {
namespace internal {

static const string token_ID = string("ID");
static const string token_CN = string("CN");
static const string token_DS = string("DS");
static const string token_DT = string("DT");
static const string token_FO = string("FO");
static const string token_KS = string("KS");
static const string token_LB = string("LB");
static const string token_PG = string("PG");
static const string token_PI = string("PI");
static const string token_PL = string("PL");
static const string token_PU = string("PU");
static const string token_SM = string("SM");

static const string feature_DQ = string("DeletionQV");
static const string feature_DT = string("DeletionTag");
static const string feature_IQ = string("InsertionQV");
static const string feature_MQ = string("MergeQV");
static const string feature_SQ = string("SubstitutionQV");
static const string feature_ST = string("SubstitutionTag");
static const string feature_IP = string("Ipd");
static const string feature_PW = string("PulseWidth");

static const string token_RT = string("READTYPE");
static const string token_BK = string("BINDINGKIT");
static const string token_SK = string("SEQUENCINGKIT");
static const string token_BV = string("BASECALLERVERSION");

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
    }
    return string();
}

static map<string, BaseFeature> nameToFeature;

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
    }
}

static inline
bool IsBaseFeature(const std::string& name)
{
    InitNameToFeature();
    return nameToFeature.find(name) != nameToFeature.cend();
}

static inline
BaseFeature FromName(const std::string& name)
{
    InitNameToFeature();
    return nameToFeature.at(name);
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

ReadGroupInfo::ReadGroupInfo(void)
    : readType_("UNKNOWN")
{ }

ReadGroupInfo::ReadGroupInfo(const std::string& id)
    : id_(id)
    , readType_("UNKNOWN")
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
    features_ = std::move(other.features_);
    return *this;
}

void ReadGroupInfo::DecodeSamDescription(const std::string& description)
{
    // split on semicolons
    // for each, split on equal
    //    determine name ->

    const vector<string>& tokens = internal::Split(description, ';');
    if (tokens.empty())
        return;

    // iterate over tokens
    auto tokenEnd  = tokens.cend();
    for (auto tokenIter = tokens.cbegin(); tokenIter != tokenEnd; ++tokenIter) {

        const string& token = *tokenIter;

        const size_t foundEqual = token.find('=');
        if (foundEqual == string::npos)
            continue;

        const string& key = token.substr(0,foundEqual);
        const string& value = token.substr(foundEqual+1);

        if      (key == internal::token_RT) readType_ = value;
        else if (key == internal::token_BK) bindingKit_ = value;
        else if (key == internal::token_BV) basecallerVersion_ = value;
        else if (key == internal::token_SK) sequencingKit_ = value;
        else if (internal::IsBaseFeature(key))
            features_[internal::FromName(key)] = value;
    }
}

std::string ReadGroupInfo::EncodeSamDescription(void) const
{
    string result;
    result.reserve(256);
    result.append(std::string("READTYPE=" + readType_));

    const auto featureEnd = features_.cend();
    auto featureIter = features_.cbegin();
    for ( ; featureIter != featureEnd; ++featureIter ) {
        const string& featureName = internal::BaseFeatureName(featureIter->first);
        if (featureName.empty() || featureIter->second.empty())
            continue;
        result.append(string(';' + featureName + '=' + featureIter->second));
    }

    if (!bindingKit_.empty())        result.append(";BINDINGKIT="+bindingKit_);
    if (!sequencingKit_.empty())     result.append(";SEQUENCINGKIT="+sequencingKit_);
    if (!basecallerVersion_.empty()) result.append(";BASECALLERVERSION="+basecallerVersion_);

    return result;
}

ReadGroupInfo ReadGroupInfo::FromSam(const string& sam)
{
    // pop off '@RG\t', then split rest of line into tokens
    const vector<string>& tokens = internal::Split(sam.substr(4), '\t');
    if (tokens.empty())
        return ReadGroupInfo();

    ReadGroupInfo rg;
    map<string, string> custom;

    for (const string& token : tokens) {
        const string& tokenTag   = token.substr(0,2);
        const string& tokenValue = token.substr(3);

        // set read group info
        if      (tokenTag == internal::token_ID) rg.Id(tokenValue);
        else if (tokenTag == internal::token_CN) rg.SequencingCenter(tokenValue);
        else if (tokenTag == internal::token_DT) rg.Date(tokenValue);
        else if (tokenTag == internal::token_FO) rg.FlowOrder(tokenValue);
        else if (tokenTag == internal::token_KS) rg.KeySequence(tokenValue);
        else if (tokenTag == internal::token_LB) rg.Library(tokenValue);
        else if (tokenTag == internal::token_PG) rg.Programs(tokenValue);
        else if (tokenTag == internal::token_PI) rg.PredictedInsertSize(tokenValue);
        else if (tokenTag == internal::token_PU) rg.MovieName(tokenValue);
        else if (tokenTag == internal::token_SM) rg.Sample(tokenValue);
        else if (tokenTag == internal::token_DS) rg.DecodeSamDescription(tokenValue);

        // otherwise, "custom" tag
        else
            custom[tokenTag] = tokenValue;
    }
    rg.CustomTags(custom);

    return rg;
}

std::string ReadGroupInfo::ToSam(void) const
{
    stringstream out;
    out << "@RG"
        << internal::MakeSamTag(internal::token_ID, id_)
        << internal::MakeSamTag(internal::token_PL, Platform());

    const string& description = EncodeSamDescription();
    if (!description.empty())
        out << internal::MakeSamTag(internal::token_DS, description);

    if (!sequencingCenter_.empty())    out << internal::MakeSamTag(internal::token_CN, sequencingCenter_);
    if (!date_.empty())                out << internal::MakeSamTag(internal::token_DT, date_);
    if (!flowOrder_.empty())           out << internal::MakeSamTag(internal::token_FO, flowOrder_);
    if (!keySequence_.empty())         out << internal::MakeSamTag(internal::token_KS, keySequence_);
    if (!library_.empty())             out << internal::MakeSamTag(internal::token_LB, library_);
    if (!programs_.empty())            out << internal::MakeSamTag(internal::token_PG, programs_);
    if (!predictedInsertSize_.empty()) out << internal::MakeSamTag(internal::token_PI, predictedInsertSize_);
    if (!movieName_.empty())           out << internal::MakeSamTag(internal::token_PU, movieName_);
    if (!sample_.empty())              out << internal::MakeSamTag(internal::token_SM, sample_);

    // append any custom tags
    map<string, string>::const_iterator customIter = custom_.cbegin();
    map<string, string>::const_iterator customEnd  = custom_.cend();
    for ( ; customIter != customEnd; ++customIter )
        out << internal::MakeSamTag(customIter->first, customIter->second);

    return out.str();
}
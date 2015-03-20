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

#ifndef READGROUPINFO_H
#define READGROUPINFO_H

#include "pbbam/Config.h"
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

enum class BaseFeature
{
    DELETION_QV
  , DELETION_TAG
  , INSERTION_QV
  , MERGE_QV
  , SUBSTITUTION_QV
  , SUBSTITUTION_TAG
  , IPD
  , PULSE_WIDTH
};

class PBBAM_EXPORT ReadGroupInfo
{
public:
    /// \name Conversion & Validation
    ///

    static ReadGroupInfo FromSam(const std::string& sam);

    static std::string ToSam(const ReadGroupInfo& rg);

    /// \}

public:
    /// \name Constructors & Related Methods
    /// \{

    ReadGroupInfo(void);
    ReadGroupInfo(const std::string& id);
    ReadGroupInfo(const ReadGroupInfo& other);
    ReadGroupInfo(ReadGroupInfo&& other);
    ReadGroupInfo& operator=(const ReadGroupInfo& other);
    ReadGroupInfo& operator=(ReadGroupInfo&& other);
    ~ReadGroupInfo(void);

    /// \}

public:
    /// \name Attributes
    /// \{

    const std::string& BasecallerVersion(void) const;

    bool HasBaseFeature(const BaseFeature& feature) const;

    std::string BaseFeatureTag(const BaseFeature& feature) const;

    std::string BindingKit(void) const;

    std::map<std::string, std::string> CustomTags(void) const;

    std::string Date(void) const;

    std::string FlowOrder(void) const;

    std::string Id(void) const;

    std::string KeySequence(void) const;

    std::string Library(void) const;

    std::string MovieName(void) const;

    std::string Platform(void) const;

    std::string PredictedInsertSize(void) const;

    std::string Programs(void) const;

    std::string ReadType(void) const;

    std::string Sample(void) const;

    std::string SequencingCenter(void) const;

    std::string SequencingKit(void) const;

    /// \}

    /// \name Conversion & Validation
    /// \{

    bool IsValid(void) const;

    std::string ToSam(void) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    ReadGroupInfo& BasecallerVersion(const std::string& versionNumber);

    ReadGroupInfo& BaseFeatureTag(const BaseFeature& feature,
                                  const std::string& tag);

    ReadGroupInfo& BindingKit(const std::string& kitNumber);

    ReadGroupInfo& CustomTags(const std::map<std::string, std::string>& custom);

    ReadGroupInfo& Date(const std::string& date);

    ReadGroupInfo& FlowOrder(const std::string& order);

    ReadGroupInfo& Id(const std::string& id);

    ReadGroupInfo& KeySequence(const std::string& sequence);

    ReadGroupInfo& Library(const std::string& library);

    ReadGroupInfo& MovieName(const std::string& id);

    ReadGroupInfo& PredictedInsertSize(const std::string& size);

    ReadGroupInfo& Programs(const std::string& programs);

    ReadGroupInfo& ReadType(const std::string& type);

    ReadGroupInfo& Sample(const std::string& sample);

    ReadGroupInfo& SequencingCenter(const std::string& center);

    ReadGroupInfo& SequencingKit(const std::string& kitNumber);

    /// \}

private:
    std::string id_;                     // ID * Unique ID required for valid SAM/BAM header *
    std::string sequencingCenter_;       // CN
    std::string date_;                   // DT * (ISO 8601) *
    std::string flowOrder_;              // FO
    std::string keySequence_;            // KS
    std::string library_;                // LB
    std::string programs_;               // PG
    std::string predictedInsertSize_;    // PI
    std::string movieName_;              // PU * more explicit, in place of "platform unit" *
    std::string sample_;                 // SM

    // DS:<Description> components
    std::string readType_;
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    std::map<BaseFeature, std::string> features_;

    // custom attributes
    std::map<std::string, std::string> custom_; // tag => value

private:
    std::string EncodeSamDescription(void) const;
    void DecodeSamDescription(const std::string& description);
};

inline const std::string& ReadGroupInfo::BasecallerVersion(void) const
{ return basecallerVersion_; }

inline ReadGroupInfo& ReadGroupInfo::BasecallerVersion(const std::string& versionNumber)
{ basecallerVersion_ = versionNumber; return *this; }

inline std::string ReadGroupInfo::BaseFeatureTag(const BaseFeature& feature) const
{
    const auto iter = features_.find(feature);
    if (iter == features_.end())
        return std::string();
    return iter->second;
}

inline ReadGroupInfo& ReadGroupInfo::BaseFeatureTag(const BaseFeature& feature,
                                                    const std::string& tag)
{ features_[feature] = tag; return *this; }

inline std::string ReadGroupInfo::BindingKit(void) const
{ return bindingKit_; }

inline ReadGroupInfo& ReadGroupInfo::BindingKit(const std::string& kitNumber)
{ bindingKit_ = kitNumber; return *this; }

inline std::map<std::string, std::string> ReadGroupInfo::CustomTags(void) const
{ return custom_; }

inline ReadGroupInfo& ReadGroupInfo::CustomTags(const std::map<std::string, std::string>& custom)
{ custom_ = custom; return *this; }

inline std::string ReadGroupInfo::Date(void) const
{ return date_; }

inline ReadGroupInfo& ReadGroupInfo::Date(const std::string& date)
{ date_ = date; return *this; }

inline std::string ReadGroupInfo::FlowOrder(void) const
{ return flowOrder_; }

inline ReadGroupInfo& ReadGroupInfo::FlowOrder(const std::string& order)
{ flowOrder_ = order; return *this; }

inline bool ReadGroupInfo::HasBaseFeature(const BaseFeature& feature) const
{ return features_.find(feature) != features_.end(); }

inline std::string ReadGroupInfo::Id(void) const
{ return id_; }

inline ReadGroupInfo& ReadGroupInfo::Id(const std::string& id)
{ id_ = id; return *this; }

inline bool ReadGroupInfo::IsValid(void) const
{ return !id_.empty(); }

inline std::string ReadGroupInfo::KeySequence(void) const
{ return keySequence_; }

inline ReadGroupInfo& ReadGroupInfo::KeySequence(const std::string& sequence)
{ keySequence_ = sequence; return *this; }

inline std::string ReadGroupInfo::Library(void) const
{ return library_; }

inline ReadGroupInfo& ReadGroupInfo::Library(const std::string& library)
{ library_ = library; return *this; }

inline std::string ReadGroupInfo::MovieName(void) const
{ return movieName_; }

inline ReadGroupInfo& ReadGroupInfo::MovieName(const std::string& movieName)
{ movieName_ = movieName; return *this; }

inline std::string ReadGroupInfo::Platform(void) const
{ return std::string("PACBIO"); }

inline std::string ReadGroupInfo::PredictedInsertSize(void) const
{ return predictedInsertSize_; }

inline ReadGroupInfo& ReadGroupInfo::PredictedInsertSize(const std::string& size)
{ predictedInsertSize_ = size; return *this; }

inline std::string ReadGroupInfo::Programs(void) const
{ return programs_; }

inline ReadGroupInfo& ReadGroupInfo::Programs(const std::string& programs)
{ programs_ = programs; return *this; }

inline std::string ReadGroupInfo::ReadType(void) const
{ return readType_; }

inline ReadGroupInfo& ReadGroupInfo::ReadType(const std::string& type)
{ readType_ = type; return *this; }

inline std::string ReadGroupInfo::Sample(void) const
{ return sample_; }

inline ReadGroupInfo& ReadGroupInfo::Sample(const std::string& sample)
{ sample_ = sample; return *this; }

inline std::string ReadGroupInfo::SequencingCenter(void) const
{ return sequencingCenter_; }

inline ReadGroupInfo& ReadGroupInfo::SequencingCenter(const std::string& center)
{ sequencingCenter_ = center; return *this; }

inline std::string ReadGroupInfo::SequencingKit(void) const
{ return sequencingKit_; }

inline ReadGroupInfo& ReadGroupInfo::SequencingKit(const std::string& kitNumber)
{ sequencingKit_ = kitNumber; return *this; }

inline std::string ReadGroupInfo::ToSam(const ReadGroupInfo& rg)
{ return rg.ToSam(); }

} // namespace BAM
} // namespace PacBio

#endif // READGROUPINFO_H

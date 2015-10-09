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

#ifndef BAMHEADER_H
#define BAMHEADER_H

#include "pbbam/Config.h"
#include "pbbam/ProgramInfo.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/SequenceInfo.h"
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { class BamHeaderPrivate; }

class PBBAM_EXPORT BamHeader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    BamHeader(void);
    BamHeader(const std::string& samHeaderText);
    BamHeader(const BamHeader& other);
    BamHeader(BamHeader&& other);
    BamHeader& operator=(const BamHeader& other);
    BamHeader& operator=(BamHeader&& other);
    ~BamHeader(void);

    BamHeader DeepCopy(void) const;

    /// \}

public:
    /// \name General
    /// \{

    std::string PacBioBamVersion(void) const;
    std::string SortOrder(void) const;
    std::string Version(void) const;

    /// \}

    /// \name Read Groups
    /// \{

    bool HasReadGroup(const std::string& id) const;
    ReadGroupInfo ReadGroup(const std::string& id) const;
    std::vector<std::string> ReadGroupIds(void) const;
    std::vector<ReadGroupInfo> ReadGroups(void) const;

    /// \}

    /// \name Sequences
    /// \{

    bool HasSequence(const std::string& name) const;
    size_t NumSequences(void) const;
    int32_t SequenceId(const std::string& name) const;
    std::string SequenceLength(const int32_t id) const;
    std::string SequenceName(const int32_t id) const;
    std::vector<std::string> SequenceNames(void) const;
    SequenceInfo Sequence(const int32_t id) const;
    SequenceInfo Sequence(const std::string& name) const;
    std::vector<SequenceInfo> Sequences(void) const;

    /// \}

    /// \name Programs
    /// \{

    bool HasProgram(const std::string& id) const;
    ProgramInfo Program(const std::string& id) const;
    std::vector<std::string> ProgramIds(void) const;
    std::vector<ProgramInfo> Programs(void) const;

    /// \}

    /// \name Comments
    /// \{

    std::vector<std::string> Comments(void) const;

    /// \}

    /// \name Conversion Methods
    /// \{

    std::string ToSam(void) const;

    /// \}

public:

    /// \name General
    /// \{

    BamHeader& PacBioBamVersion(const std::string& version);
    BamHeader& SortOrder(const std::string& order);
    BamHeader& Version(const std::string& version);

    /// \}

    /// \name Read Groups
    /// \{

    BamHeader& AddReadGroup(const ReadGroupInfo& readGroup);
    BamHeader& ClearReadGroups(void);
    BamHeader& ReadGroups(const std::vector<ReadGroupInfo>& readGroups);

    /// \}

    /// \name Sequences
    /// \{

    BamHeader& AddSequence(const SequenceInfo& sequence);
    BamHeader& ClearSequences(void);
    BamHeader& Sequences(const std::vector<SequenceInfo>& sequences);

    /// \}

    /// \name Programs
    /// \{

    BamHeader& AddProgram(const ProgramInfo& pg);
    BamHeader& ClearPrograms(void);
    BamHeader& Programs(const std::vector<ProgramInfo>& programs);

    /// \}

    /// \name Comments
    /// \{

    BamHeader& AddComment(const std::string& comment);
    BamHeader& ClearComments(void);
    BamHeader& Comments(const std::vector<std::string>& comments);

    /// \}

private:
    PBBAM_SHARED_PTR<internal::BamHeaderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#endif // BAMHEADER_H

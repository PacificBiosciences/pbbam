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
/// \file BamHeader.h
/// \brief Defines the BamHeader class.
//
// Author: Derek Barnett

#ifndef BAMHEADER_H
#define BAMHEADER_H

#include "pbbam/Config.h"
#include "pbbam/ProgramInfo.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/SequenceInfo.h"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { class BamHeaderPrivate; }

/// \brief The BamHeader class represents the header section of the %BAM file.
///
/// It provides metadata about the file including file version, reference
/// sequences, read groups, comments, etc.
///
/// A BamHeader may be fetched from a BamFile to view an existing file's
/// metadata. Or one may be created/edited for use with writing to a new file
/// (via BamWriter).
///
/// \note A particular BamHeader is likely to be re-used in lots of places
///       throughout the library, for read-only purposes. For this reason, even
///       though a BamHeader may be returned by value, it is essentially a thin
///       wrapper for a shared-pointer to the actual data. This means, though,
///       that if you need to edit an existing BamHeader for use with a
///       BamWriter, please consider using BamHeader::DeepCopy. Otherwise any
///       modifications will affect all BamHeaders that are sharing its
///       underlying data.
///
class PBBAM_EXPORT BamHeader
{
public:
    /// \name Constructors & Related Methods
    /// \{

    ///
    /// \brief Creates an empty BamHeader
    ///
    BamHeader();

    ///
    /// \brief Creates a BamHeader from SAM-formatted text
    /// \param samHeaderText
    ///
    BamHeader(const std::string& samHeaderText);

    BamHeader(const BamHeader& other) = default;
    BamHeader(BamHeader&& other) = default;
    BamHeader& operator=(const BamHeader& other) = default;
    BamHeader& operator=(BamHeader&& other) = default;
    ~BamHeader() = default;

    /// \brief Detaches underlying data from the shared-pointer, returning a
    ///        independent copy of the header contents.
    ///
    /// This ensures that any modifications to the newly returned BamHeader do
    /// not affect other BamHeader objects that were sharing its underlying data.
    ///
    BamHeader DeepCopy() const;

    /// \}

public:
    /// \name Operators
    /// \{

    /// \brief Merges another header with this one.
    ///
    /// Headers must be compatible for merging. This means that their Version,
    /// SortOrder, PacBioBamVersion (and in the case of aligned BAM data,
    /// Sequences) must all match. If not, an exception will be thrown.
    ///
    /// \param[in] other  header to merge with this one
    /// \returns reference to this header
    ///
    /// \throws std::runtime_error if the headers are not compatible
    ///
    BamHeader& operator+=(const BamHeader& other);

    /// \brief Creates a new, merged header.
    ///
    /// Headers must be compatible for merging. This means that their Version,
    /// SortOrder, PacBioBamVersion (and in the case of aligned BAM data,
    /// Sequences) must all match. If not, an exception will be thrown.
    ///
    /// Both original headers (this header and \p other) will not be modified.
    ///
    /// \param[in] other  header to merge with this one
    /// \returns merged header
    ///
    /// \throws std::runtime_error if the headers are not compatible
    ///
    BamHeader operator+(const BamHeader& other) const;

    /// \}

public:
    /// \name General Attributes
    /// \{

    /// \returns the %PacBio %BAM version number (\@HD:pb)
    ///
    /// \note This is different from the SAM/BAM version number
    /// \sa BamHeader::Version.
    ///
    std::string PacBioBamVersion() const;

    /// \returns the sort order used
    ///
    /// Valid values: "unknown", "unsorted", "queryname", or "coordinate"
    ///
    std::string SortOrder() const;

    /// \returns the SAM/BAM version number (\@HD:VN)
    ///
    /// \note This is different from the %PacBio %BAM version number
    /// \sa BamHeader::PacBioBamVersion
    ///
    std::string Version() const;

    /// \}

public:
    /// \name Read Groups
    /// \{

    /// \returns true if the header contains a read group with \p id (\@RG:ID)
    bool HasReadGroup(const std::string& id) const;

    /// \returns a ReadGroupInfo object representing the read group matching
    ///          \p id (\@RG:ID)
    /// \throws std::runtime_error if \p id is unknown
    ///
    ReadGroupInfo ReadGroup(const std::string& id) const;

    /// \returns vector of read group IDs listed in this header
    std::vector<std::string> ReadGroupIds() const;

    /// \returns vector of ReadGroupInfo objects, representing all read groups
    ///          listed in this header
    ///
    std::vector<ReadGroupInfo> ReadGroups() const;

    /// \}

public:
    /// \name Sequences
    /// \{

    /// \returns true if header contains a sequence with \p name (\@SQ:SN)
    bool HasSequence(const std::string& name) const;

    /// \returns number of sequences (\@SQ entries) stored in this header
    size_t NumSequences() const;

    /// \returns numeric ID for sequence matching \p name (\@SQ:SN)
    ///
    /// This is the numeric ID used elsewhere throughout the API.
    ///
    /// \throws std::runtime_error if \p name is unknown
    /// \sa BamReader::ReferenceId, PbiReferenceIdFilter,
    ///     PbiRawMappedData::tId_
    ///
    int32_t SequenceId(const std::string& name) const;

    /// \returns the length of the sequence (\@SQ:LN, e.g. chromosome length) at
    ///          index \p id
    ///
    /// \sa SequenceInfo::Length, BamHeader::SequenceId
    ///
    std::string SequenceLength(const int32_t id) const;

    /// \returns the name of the sequence (\@SQ:SN) at index \p id
    ///
    /// \sa SequenceInfo::Name, BamHeader::SequenceId
    ///
    std::string SequenceName(const int32_t id) const;

    /// \returns vector of sequence names (\@SQ:SN) stored in this header
    ///
    /// Position in the vector is equivalent to SequenceId.
    ///
    std::vector<std::string> SequenceNames() const;

    /// \returns SequenceInfo object at index \p id
    ///
    /// \throws std::out_of_range if \p is an invalid or unknown index
    /// \sa BamHeader::SequenceId
    ///
    SequenceInfo Sequence(const int32_t id) const;

    /// \returns SequenceInfo for the sequence matching \p name
    SequenceInfo Sequence(const std::string& name) const;

    /// \returns vector of SequenceInfo objects representing the sequences
    ///          (\@SQ entries) stored in this header
    ///
    std::vector<SequenceInfo> Sequences() const;

    /// \}

public:
    /// \name Programs
    /// \{

    /// \returns true if this header contains a program entry with ID (\@PG:ID)
    ///          matching \p id
    ///
    bool HasProgram(const std::string& id) const;

    /// \returns ProgramInfo object for the program entry matching \p id
    /// \throws std::runtime_error if \p id is unknown
    ///
    ProgramInfo Program(const std::string& id) const;

    /// \returns vector of program IDs (\@PG:ID)
    std::vector<std::string> ProgramIds() const;

    /// \returns vector of ProgramInfo objects representing program entries
    ///          (\@PG) stored in this heder
    ///
    std::vector<ProgramInfo> Programs() const;

    /// \}

public:
    /// \name Comments
    /// \{

    /// \returns vector of comment (\@CO) strings
    std::vector<std::string> Comments() const;

    /// \}

public:
    /// \name Conversion Methods
    /// \{

    /// \returns SAM-header-formatted string representing this header's data
    std::string ToSam() const;

    /// \}

public:

    /// \name General Attributes
    /// \{

    /// \brief Sets this header's PacBioBAM version number (\@HD:pb).
    ///
    /// \returns reference to this object
    /// \throws std::runtime_error if version number cannot be parsed or
    ///         is less than the minimum version allowed.
    ///
    BamHeader& PacBioBamVersion(const std::string& version);

    /// \brief Sets this header's sort order label (\@HD:SO).
    ///
    /// Valid values: "unknown", "unsorted", "queryname", or "coordinate"
    ///
    /// \returns reference to this object
    ///
    BamHeader& SortOrder(const std::string& order);

    /// \brief Sets this header's SAM/BAM version number (\@HD:VN).
    ///
    /// \returns reference to this object
    ///
    BamHeader& Version(const std::string& version);

    /// \}

public:
    /// \name Read Groups
    /// \{

    /// \brief Appends a read group entry (\@RG) to this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& AddReadGroup(const ReadGroupInfo& readGroup);

    /// \brief Removes all read group entries from this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& ClearReadGroups();

    /// \brief Replaces this header's list of read group entries with those in
    ///        \p readGroups.
    ///
    /// \returns reference to this object
    ///
    BamHeader& ReadGroups(const std::vector<ReadGroupInfo>& readGroups);

    /// \}

public:
    /// \name Sequences
    /// \{

    /// \brief Appends a sequence entry (\@SQ) to this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& AddSequence(const SequenceInfo& sequence);

    /// \brief Removes all sequence entries from this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& ClearSequences();

    /// \brief Replaces this header's list of sequence entries with those in
    ///       \p sequences.
    ///
    /// \returns reference to this object
    ///
    BamHeader& Sequences(const std::vector<SequenceInfo>& sequences);

    /// \}

public:
    /// \name Programs
    /// \{

    /// \brief Appends a program entry (\@PG) to this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& AddProgram(const ProgramInfo& pg);

    /// \brief Removes all program entries from this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& ClearPrograms();

    /// \brief Replaces this header's list of program entries with those in
    ///        \p programs.
    ///
    /// \returns reference to this object
    ///
    BamHeader& Programs(const std::vector<ProgramInfo>& programs);

    /// \}

public:
    /// \name Comments
    /// \{

    /// \brief Appends a comment (\@CO) to this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& AddComment(const std::string& comment);

    /// \brief Removes all comments from this header.
    ///
    /// \returns reference to this object
    ///
    BamHeader& ClearComments();

    /// \brief Replaces this header's list of comments with those in \p comments.
    ///
    /// \returns reference to this object
    ///
    BamHeader& Comments(const std::vector<std::string>& comments);

    /// \}

private:
    std::shared_ptr<internal::BamHeaderPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#include "pbbam/internal/BamHeader.inl"

#endif // BAMHEADER_H

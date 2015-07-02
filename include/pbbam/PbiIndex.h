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
// Author: Derek Barnett

#ifndef PBIINDEX_H
#define PBIINDEX_H

#include "pbbam/Config.h"
#include "pbbam/LocalContextFlags.h"
#include "pbbam/PbiFile.h"
#include "pbbam/Strand.h"
#include <deque>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal { class PbiIndexPrivate; }

enum class SubreadField
{
    RG_ID
  , Q_START
  , Q_END
  , ZMW
  , READ_QUALITY
  , VIRTUAL_OFFSET
};

enum class MappedField
{
    T_ID
  , T_START
  , T_END
  , A_START
  , A_END
  , N_M
  , N_MM
  , N_INS
  , N_DEL
  , MAP_QUALITY
  , STRAND
};

enum class BarcodeField
{
    BC_LEFT
  , BC_RIGHT
  , BC_QUALITY
  , CONTEXT_FLAG
};

enum class CompareType
{
    EQUAL
  , LESS_THAN
  , LESS_THAN_EQUAL
  , GREATER_THAN
  , GREATER_THAN_EQUAL
  , NOT_EQUAL
};

//
// Contiguous reads that satisfy a query will be returned as a block.
// This is to help minimize number of seeks (or even unneccesary checks).
//
// An index query can iterate over the lookup result 'IndexResultBlocks' list to
// perform a seek and fetch 'numReads' consecutive records before needing to
// seek again.
//
struct PBBAM_EXPORT IndexResultBlock
{
public:
    IndexResultBlock(void);
    IndexResultBlock(size_t idx, size_t numReads);

public:
    bool operator==(const IndexResultBlock& other) const;
    bool operator!=(const IndexResultBlock& other) const;

public:
    size_t  firstIndex_;
    size_t  numReads_;
    int64_t virtualOffset_;
};

typedef std::deque<IndexResultBlock> IndexResultBlocks;

typedef std::vector<size_t>       IndexList;
typedef std::pair<size_t, size_t> IndexRange;

template<typename FieldType, typename ValueType>
struct IndexRequestBase
{
public:
    FieldType field_;
    ValueType value_;
    CompareType compareType_;

protected:
    IndexRequestBase(const FieldType field,
                     const ValueType& value,
                     const CompareType compareType = CompareType::EQUAL);
};

// all multi-requests use CompareType::EQUAL
template<typename FieldType, typename ValueType>
struct IndexMultiRequestBase
{
public:
    FieldType field_;
    std::vector<ValueType> values_;

protected:
    IndexMultiRequestBase(const FieldType field,
                          const std::vector<ValueType>& values);
};

class PBBAM_EXPORT PbiIndex
{
public:
    /// \name Constructors & Related Methods
    /// \{

    PbiIndex(const std::string& pbiFilename);
    PbiIndex(const PbiIndex& other);
    PbiIndex(PbiIndex&& other);
    PbiIndex& operator=(const PbiIndex& other);
    PbiIndex& operator=(PbiIndex&& other);
    ~PbiIndex(void);

    /// \}

public:
    // PBI attributes
    bool HasBarcodeData(void) const;
    bool HasMappedData(void) const;
    bool HasReferenceData(void) const;
    bool HasSection(const PbiFile::Section section) const;

    PbiFile::Sections FileSections(void) const;
    uint32_t NumReads(void) const;
    PbiFile::VersionEnum Version(void) const;

public:

    template<typename FieldType, typename ValueType>
    IndexList RawIndices(const IndexRequestBase<FieldType, ValueType>& request) const;

    template<typename FieldType, typename ValueType>
    IndexList RawIndices(const IndexMultiRequestBase<FieldType, ValueType>& request) const;

    template<typename FieldType, typename ValueType>
    IndexResultBlocks Lookup(const IndexRequestBase<FieldType, ValueType>& request) const;

    template<typename FieldType, typename ValueType>
    IndexResultBlocks Lookup(const IndexMultiRequestBase<FieldType, ValueType>& request) const;

    IndexResultBlocks LookupReference(const int32_t tId) const;

    const std::vector<int64_t>& VirtualFileOffsets(void) const;

private:
    PbiIndex(void);
    std::unique_ptr<internal::PbiIndexPrivate> d_;
};

template<SubreadField field, typename ValueType>
class SubreadIndexRequest : public IndexRequestBase<SubreadField, ValueType>
{
public:
    SubreadIndexRequest(const ValueType& value,
                        const CompareType& compareType = CompareType::EQUAL);
};

template<SubreadField field, typename ValueType>
class SubreadIndexMultiRequest : public IndexMultiRequestBase<SubreadField, ValueType>
{
public:
    SubreadIndexMultiRequest(const std::vector<ValueType>& values);
};

typedef SubreadIndexRequest<SubreadField::RG_ID,        int32_t>  ReadGroupIndexRequest;
typedef SubreadIndexRequest<SubreadField::Q_START,      int32_t>  QueryStartIndexRequest;
typedef SubreadIndexRequest<SubreadField::Q_END,        int32_t>  QueryEndIndexRequest;
typedef SubreadIndexRequest<SubreadField::ZMW,          int32_t>  ZmwIndexRequest;
typedef SubreadIndexRequest<SubreadField::READ_QUALITY, uint16_t> ReadQualityIndexRequest;

typedef SubreadIndexMultiRequest<SubreadField::RG_ID,        int32_t>  ReadGroupIndexMultiRequest;
typedef SubreadIndexMultiRequest<SubreadField::Q_START,      int32_t>  QueryStartIndexMultiRequest;
typedef SubreadIndexMultiRequest<SubreadField::Q_END,        int32_t>  QueryEndIndexMultiRequest;
typedef SubreadIndexMultiRequest<SubreadField::ZMW,          int32_t>  ZmwIndexMultiRequest;
typedef SubreadIndexMultiRequest<SubreadField::READ_QUALITY, uint16_t> ReadQualityIndexMultiRequest;

template<MappedField field, typename ValueType>
class MappedIndexRequest : public IndexRequestBase<MappedField, ValueType>
{
public:
    MappedIndexRequest(const ValueType& value, const
                       CompareType& compareType = CompareType::EQUAL);
};

template<MappedField field, typename ValueType>
class MappedIndexMultiRequest : public IndexMultiRequestBase<SubreadField, ValueType>
{
public:
    MappedIndexMultiRequest(const std::vector<ValueType>& values);
};

typedef MappedIndexRequest<MappedField::T_ID,         int32_t> ReferenceIdIndexRequest;
typedef MappedIndexRequest<MappedField::T_START,      int32_t> ReferenceStartIndexRequest;
typedef MappedIndexRequest<MappedField::T_END,        int32_t> ReferenceEndIndexRequest;
typedef MappedIndexRequest<MappedField::A_START,      int32_t> AlignedStartIndexRequest;
typedef MappedIndexRequest<MappedField::A_END,        int32_t> AlignedEndIndexRequest;
typedef MappedIndexRequest<MappedField::N_M,          int32_t> NumMatchesIndexRequest;
typedef MappedIndexRequest<MappedField::N_MM,         int32_t> NumMismatchesIndexRequest;
typedef MappedIndexRequest<MappedField::N_INS,        int32_t> NumInsertionsIndexRequest;
typedef MappedIndexRequest<MappedField::N_DEL,        int32_t> NumDeletionsIndexRequest;
typedef MappedIndexRequest<MappedField::MAP_QUALITY,  uint8_t> MapQualityIndexRequest;
typedef MappedIndexRequest<MappedField::STRAND,       Strand>  StrandIndexRequest;

typedef MappedIndexMultiRequest<MappedField::T_ID,         int32_t> ReferenceIdIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::T_START,      int32_t> ReferenceStartIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::T_END,        int32_t> ReferenceEndIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::A_START,      int32_t> AlignedStartIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::A_END,        int32_t> AlignedEndIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::N_M,          int32_t> NumMatchesIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::N_MM,         int32_t> NumMismatchesIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::N_INS,        int32_t> NumInsertionsIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::N_DEL,        int32_t> NumDeletionsIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::MAP_QUALITY,  uint8_t> MapQualityIndexMultiRequest;
typedef MappedIndexMultiRequest<MappedField::STRAND,       Strand>  StrandIndexMultiRequest;

template<BarcodeField field, typename ValueType>
class BarcodeIndexRequest : public IndexRequestBase<BarcodeField, ValueType>
{
public:
    BarcodeIndexRequest(const ValueType& value,
                        const CompareType& compareType = CompareType::EQUAL);
};

template<BarcodeField field, typename ValueType>
class BarcodeIndexMultiRequest : public IndexMultiRequestBase<BarcodeField, ValueType>
{
public:
    BarcodeIndexMultiRequest(const std::vector<ValueType>& values);
};

typedef BarcodeIndexRequest<BarcodeField::BC_LEFT,      uint16_t> BarcodeLeftIndexRequest;
typedef BarcodeIndexRequest<BarcodeField::BC_RIGHT,     uint16_t> BarcodeRightIndexRequest;
typedef BarcodeIndexRequest<BarcodeField::BC_QUALITY,   uint8_t>  BarcodeQualityIndexRequest;
typedef BarcodeIndexRequest<BarcodeField::CONTEXT_FLAG, LocalContextFlags> ContextFlagIndexRequest;

typedef BarcodeIndexMultiRequest<BarcodeField::BC_LEFT,      uint16_t> BarcodeLeftIndexMultiRequest;
typedef BarcodeIndexMultiRequest<BarcodeField::BC_RIGHT,     uint16_t> BarcodeRightIndexMultiRequest;
typedef BarcodeIndexMultiRequest<BarcodeField::BC_QUALITY,   uint8_t>  BarcodeQualityIndexMultiRequest;
typedef BarcodeIndexMultiRequest<BarcodeField::CONTEXT_FLAG, LocalContextFlags> ContextFlagIndexMultiRequest;

} // namespace BAM
} // namespace PacBio

#include "internal/PbiIndex_p.inl"

#endif // PBIINDEX_H

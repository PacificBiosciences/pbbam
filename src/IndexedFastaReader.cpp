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
/// \file IndexedFastaReader.cpp
/// \brief Implements the IndexedFastaReader class.
//
// Author: David Alexander

#include "pbbam/IndexedFastaReader.h"

#include "pbbam/BamRecord.h"
#include "pbbam/GenomicInterval.h"
#include "pbbam/Orientation.h"
#include "SequenceUtils.h"
#include <htslib/faidx.h>
#include <iostream>
#include <cstdlib>

namespace PacBio {
namespace BAM {

IndexedFastaReader::IndexedFastaReader(const std::string& filename)
{
    Open(filename);
}

IndexedFastaReader::IndexedFastaReader(const IndexedFastaReader& src)
{
    if (!Open(src.filename_))
        throw std::runtime_error("Cannot open file " + src.filename_);
}

IndexedFastaReader& IndexedFastaReader::operator=(const IndexedFastaReader& rhs)
{
    if (&rhs == this)
        return *this;

    Open(rhs.filename_);
    return *this;
}

IndexedFastaReader::~IndexedFastaReader(void)
{
    Close();
}

bool IndexedFastaReader::Open(const std::string &filename)
{
    faidx_t* handle = fai_load(filename.c_str());
    if (handle == nullptr)
        return false;
    else
    {
        filename_ = filename;
        handle_ = handle;
        return true;
    }
}

void IndexedFastaReader::Close(void)
{
    filename_ = "";
    if (handle_ != nullptr)
        fai_destroy(handle_);
    handle_ = nullptr;
}

#define REQUIRE_FAIDX_LOADED if (handle_ == nullptr) throw std::exception()

std::string IndexedFastaReader::Subsequence(const std::string& id,
                                            Position begin,
                                            Position end) const
{
    REQUIRE_FAIDX_LOADED;

    int len;
    // Derek: *Annoyingly* htslib seems to interpret "end" as inclusive in
    // faidx_fetch_seq, whereas it considers it exclusive in the region spec in
    // fai_fetch.  Can you please verify?
    char* rawSeq = faidx_fetch_seq(handle_, id.c_str(), begin, end - 1, &len);
    if (rawSeq == nullptr)
        throw std::runtime_error("could not fetch FASTA sequence");
    else {
        std::string seq(rawSeq);
        free(rawSeq);
        return seq;
    }
}

std::string IndexedFastaReader::Subsequence(const GenomicInterval& interval) const
{
    REQUIRE_FAIDX_LOADED;
    return Subsequence(interval.Name(), interval.Start(), interval.Stop());
}

std::string IndexedFastaReader::Subsequence(const char *htslibRegion) const
{
    REQUIRE_FAIDX_LOADED;

    int len;
    char* rawSeq = fai_fetch(handle_, htslibRegion, &len);
    if (rawSeq == nullptr)
        throw std::runtime_error("could not fetch FASTA sequence");
    else {
        std::string seq(rawSeq);
        free(rawSeq);
        return seq;
    }
}


std::string
IndexedFastaReader::ReferenceSubsequence(const BamRecord& bamRecord,
                                         const Orientation orientation,
                                         const bool gapped,
                                         const bool exciseSoftClips) const
{
    REQUIRE_FAIDX_LOADED;

    std::string subseq = Subsequence(bamRecord.ReferenceName(),
                                     bamRecord.ReferenceStart(),
                                     bamRecord.ReferenceEnd());
    const auto reverse = orientation != Orientation::GENOMIC &&
                         bamRecord.Impl().IsReverseStrand();

    if (bamRecord.Impl().IsMapped() && gapped)
    {
        size_t seqIndex = 0;
        const Cigar& cigar = bamRecord.Impl().CigarData();
        Cigar::const_iterator cigarIter = cigar.cbegin();
        Cigar::const_iterator cigarEnd = cigar.cend();
        for (; cigarIter != cigarEnd; ++cigarIter)
        {
            const CigarOperation& op = (*cigarIter);
            const CigarOperationType& type = op.Type();

            // do nothing for hard clips
            if (type != CigarOperationType::HARD_CLIP)
            {
                const size_t opLength = op.Length();

                // maybe remove soft clips
                if (type == CigarOperationType::SOFT_CLIP)
                {
                    if (!exciseSoftClips)
                    {
                        subseq.reserve(subseq.size() + opLength);
                        subseq.insert(seqIndex, opLength, '-');
                        seqIndex += opLength;
                    }
                }

                // for non-clipping operations
                else {

                    // maybe add gaps/padding
                    if (type == CigarOperationType::INSERTION)
                    {
                        subseq.reserve(subseq.size() + opLength);
                        subseq.insert(seqIndex, opLength, '-');
                    }
                    else if (type == CigarOperationType::PADDING)
                    {
                        subseq.reserve(subseq.size() + opLength);
                        subseq.insert(seqIndex, opLength, '*');
                    }

                    // update index
                    seqIndex += opLength;
                }
            }
        }
    }

    if (reverse)
        internal::ReverseComplementCaseSens(subseq);

    return subseq;
}


int IndexedFastaReader::NumSequences(void) const
{
    REQUIRE_FAIDX_LOADED;
    return faidx_nseq(handle_);
}

bool IndexedFastaReader::HasSequence(const std::string& name) const
{
    REQUIRE_FAIDX_LOADED;
    return (faidx_has_seq(handle_, name.c_str()) != 0);
}

int IndexedFastaReader::SequenceLength(const std::string& name) const
{
    REQUIRE_FAIDX_LOADED;
    int len = faidx_seq_len(handle_, name.c_str());
    if (len < 0)
        throw std::runtime_error("could not determine FASTA sequence length");
    else return len;
}

}}  // PacBio::BAM

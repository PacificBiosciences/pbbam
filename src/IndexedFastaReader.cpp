// File Description
/// \file IndexedFastaReader.cpp
/// \brief Implements the IndexedFastaReader class.
//
// Author: David Alexander

#include "PbbamInternalConfig.h"

#include "pbbam/IndexedFastaReader.h"

#include <memory>
#include <sstream>
#include <stdexcept>

#include <htslib/faidx.h>

#include "pbbam/BamRecord.h"
#include "pbbam/GenomicInterval.h"
#include "pbbam/Orientation.h"
#include "pbbam/StringUtilities.h"

#include "MemoryUtils.h"
#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {

namespace {

void ClipAndGapify(std::string& subseq, const Cigar& cigar, bool exciseSoftClips)
{
    size_t seqIndex = 0;
    for (const auto& op : cigar) {
        const auto type = op.Type();
        const auto opLength = op.Length();

        // do nothing for hard clips
        if (type == CigarOperationType::HARD_CLIP) continue;

        // maybe remove soft clips
        if (type == CigarOperationType::SOFT_CLIP) {
            if (!exciseSoftClips) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '-');
                seqIndex += opLength;
            }
        }

        // for non-clipping operations
        else {
            // maybe add gaps/padding
            if (type == CigarOperationType::INSERTION) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '-');
            } else if (type == CigarOperationType::PADDING) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '*');
            }

            // update index
            seqIndex += opLength;
        }
    }
}

}  // namespace

class IndexedFastaReader::IndexedFastaReaderPrivate
{
public:
    IndexedFastaReaderPrivate(std::string filename)
        : fastaFilename_{std::move(filename)}, faiFilename_{fastaFilename_ + ".fai"}
    {
        handle_.reset(fai_load(fastaFilename_.c_str()));
        if (handle_ == nullptr) {
            throw std::runtime_error{
                "IndexedFastaReader: could not open index file (*.fai) for FASTA file: " +
                fastaFilename_};
        }
    }

    std::string fastaFilename_;
    std::string faiFilename_;
    std::unique_ptr<faidx_t, HtslibFastaIndexDeleter> handle_;
};

IndexedFastaReader::IndexedFastaReader(std::string filename)
    : d_{std::make_unique<IndexedFastaReaderPrivate>(std::move(filename))}
{
}

IndexedFastaReader::IndexedFastaReader(const IndexedFastaReader& other)
    : d_{std::make_unique<IndexedFastaReaderPrivate>(other.d_->fastaFilename_)}
{
}

IndexedFastaReader::IndexedFastaReader(IndexedFastaReader&&) noexcept = default;

IndexedFastaReader& IndexedFastaReader::operator=(const IndexedFastaReader& rhs)
{
    IndexedFastaReader copy{rhs};
    *this = std::move(copy);
    return *this;
}

IndexedFastaReader& IndexedFastaReader::operator=(IndexedFastaReader&&) noexcept = default;

IndexedFastaReader::~IndexedFastaReader() = default;

std::string IndexedFastaReader::Subsequence(const std::string& id, Position begin,
                                            Position end) const
{
    assert(begin <= end);
    // htslib is dumb and will not consider empty intervals valid,
    // that is, a call to faidx_fetch_seq will *always* return a
    // sequence consisting of at least one base.
    if (begin == end) return std::string{};

    int len;
    // Derek: *Annoyingly* htslib seems to interpret "end" as inclusive in
    // faidx_fetch_seq, whereas it considers it exclusive in the region spec in
    // fai_fetch.  Can you please verify?
    const std::unique_ptr<char, FreeDeleter> rawSeq{
        faidx_fetch_seq(d_->handle_.get(), id.c_str(), begin, end - 1, &len)};
    if (rawSeq == nullptr) {
        std::ostringstream s;
        s << "IndexedFastaReader: could not fetch sequence from region: " << id << " [" << begin
          << ", " << end << ") in FASTA file: " << d_->fastaFilename_;
        throw std::runtime_error{s.str()};
    }
    return RemoveAllWhitespace(rawSeq.get());
}

std::string IndexedFastaReader::Subsequence(const GenomicInterval& interval) const
{
    return Subsequence(interval.Name(), interval.Start(), interval.Stop());
}

std::string IndexedFastaReader::Subsequence(const char* htslibRegion) const
{
    int len = 0;
    const std::unique_ptr<char, FreeDeleter> rawSeq(
        fai_fetch(d_->handle_.get(), htslibRegion, &len));
    if (rawSeq == nullptr) {
        throw std::runtime_error{"IndexedFastaReader: could not fetch sequence from region: " +
                                 std::string{htslibRegion} + " in FASTA file: " +
                                 d_->fastaFilename_};
    }
    return RemoveAllWhitespace(rawSeq.get());
}

std::string IndexedFastaReader::ReferenceSubsequence(const BamRecord& bamRecord,
                                                     const Orientation orientation,
                                                     const bool gapped,
                                                     const bool exciseSoftClips) const
{
    std::string subseq = Subsequence(bamRecord.ReferenceName(), bamRecord.ReferenceStart(),
                                     bamRecord.ReferenceEnd());

    if (bamRecord.Impl().IsMapped() && gapped)
        ClipAndGapify(subseq, bamRecord.Impl().CigarData(), exciseSoftClips);

    const auto reverse =
        (orientation != Orientation::GENOMIC) && bamRecord.Impl().IsReverseStrand();
    if (reverse) ReverseComplementCaseSens(subseq);

    return subseq;
}

int IndexedFastaReader::NumSequences() const { return faidx_nseq(d_->handle_.get()); }

std::vector<std::string> IndexedFastaReader::Names() const
{
    std::vector<std::string> names;
    names.reserve(NumSequences());
    for (int i = 0; i < NumSequences(); ++i)
        names.emplace_back(faidx_iseq(d_->handle_.get(), i));
    return names;
}

std::string IndexedFastaReader::Name(const size_t idx) const
{
    if (static_cast<int>(idx) >= NumSequences()) {
        std::ostringstream s;
        s << "IndexedFastaReader: cannot fetch sequence name. Index (" << idx
          << ") is larger than the number of sequences: (" << NumSequences()
          << ") in FASTA file: " << d_->fastaFilename_;
        throw std::runtime_error{s.str()};
    }
    return {faidx_iseq(d_->handle_.get(), idx)};
}

bool IndexedFastaReader::HasSequence(const std::string& name) const
{
    return (faidx_has_seq(d_->handle_.get(), name.c_str()) != 0);
}

int IndexedFastaReader::SequenceLength(const std::string& name) const
{
    const auto len = faidx_seq_len(d_->handle_.get(), name.c_str());
    if (len < 0) {
        throw std::runtime_error{"IndexedFastaReader: could not determine sequence length of " +
                                 name + " in FASTA file: " + d_->fastaFilename_};
    }
    return len;
}
}
}  // PacBio::BAM

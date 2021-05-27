#include "PbbamInternalConfig.h"

#include <pbbam/IndexedFastaReader.h>

#include <cassert>

#include <memory>
#include <sstream>
#include <stdexcept>

#include <htslib/faidx.h>
#include <pbcopper/utility/Deleters.h>

#include <pbbam/BamRecord.h>
#include <pbbam/Deleters.h>
#include <pbbam/GenomicInterval.h>
#include <pbbam/Orientation.h>
#include <pbbam/StringUtilities.h>

#include "ErrnoReason.h"
#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {

namespace {

void ClipAndGapify(std::string& subseq, const Data::Cigar& cigar, bool exciseSoftClips)
{
    size_t seqIndex = 0;
    for (const auto& op : cigar) {
        const auto type = op.Type();
        const auto opLength = op.Length();

        // do nothing for hard clips
        if (type == Data::CigarOperationType::HARD_CLIP) continue;

        // maybe remove soft clips
        if (type == Data::CigarOperationType::SOFT_CLIP) {
            if (!exciseSoftClips) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '-');
                seqIndex += opLength;
            }
        }

        // for non-clipping operations
        else {
            // maybe add gaps/padding
            if (type == Data::CigarOperationType::INSERTION) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '-');
            } else if (type == Data::CigarOperationType::PADDING) {
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
        : fastaFilename_{std::move(filename)}
        , faiFilename_{fastaFilename_ + ".fai"}
        , handle_{fai_load3(fastaFilename_.c_str(), faiFilename_.c_str(), NULL,
                            0  // do not create FAI automagically
                            )}
    {
        if (!handle_) {
            std::ostringstream msg;
            msg << "[pbbam] FASTA reader ERROR: could not load FAI index data:\n"
                << "  FASTA file: " << fastaFilename_ << '\n'
                << "  FAI file: " << faiFilename_;
            MaybePrintErrnoReason(msg);
            throw std::runtime_error{msg.str()};
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
    if (this != &rhs) *this = IndexedFastaReader{rhs};
    return *this;
}

IndexedFastaReader& IndexedFastaReader::operator=(IndexedFastaReader&&) noexcept = default;

IndexedFastaReader::~IndexedFastaReader() = default;

std::string IndexedFastaReader::Subsequence(const std::string& id, Data::Position begin,
                                            Data::Position end) const
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
    const std::unique_ptr<char, Utility::FreeDeleter> rawSeq{
        faidx_fetch_seq(d_->handle_.get(), id.c_str(), begin, end - 1, &len)};
    if (rawSeq == nullptr) {
        std::ostringstream s;
        s << "[pbbam] indexed FASTA reader ERROR: could not fetch subsequence from region: " << id
          << " [" << begin << ", " << end << ")\n"
          << "  FASTA file: " << d_->fastaFilename_;
        throw std::runtime_error{s.str()};
    }
    return RemoveAllWhitespace(rawSeq.get());
}

std::string IndexedFastaReader::Subsequence(const Data::GenomicInterval& interval) const
{
    return Subsequence(interval.Name(), interval.Start(), interval.Stop());
}

std::string IndexedFastaReader::Subsequence(const char* htslibRegion) const
{
    int len = 0;
    const std::unique_ptr<char, Utility::FreeDeleter> rawSeq(
        fai_fetch(d_->handle_.get(), htslibRegion, &len));
    if (rawSeq == nullptr) {
        std::ostringstream s;
        s << "[pbbam] indexed FASTA reader ERROR: could not fetch subsequence from region: "
          << htslibRegion << '\n'
          << "  FASTA file: " << d_->fastaFilename_;
        throw std::runtime_error{s.str()};
    }
    return RemoveAllWhitespace(rawSeq.get());
}

std::string IndexedFastaReader::ReferenceSubsequence(const BamRecord& bamRecord,
                                                     const Data::Orientation orientation,
                                                     const bool gapped,
                                                     const bool exciseSoftClips) const
{
    std::string subseq = Subsequence(bamRecord.ReferenceName(), bamRecord.ReferenceStart(),
                                     bamRecord.ReferenceEnd());

    if (bamRecord.Impl().IsMapped() && gapped)
        ClipAndGapify(subseq, bamRecord.Impl().CigarData(), exciseSoftClips);

    const auto reverse =
        (orientation != Data::Orientation::GENOMIC) && bamRecord.Impl().IsReverseStrand();
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
        s << "[pbbam] indexed FASTA reader ERROR: cannot fetch sequence name. Index (" << idx
          << ") is larger than the number of sequences: (" << NumSequences() << ")\n"
          << "  FASTA file: " << d_->fastaFilename_;
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
        std::ostringstream s;
        s << "[pbbam] indexed FASTA reader ERROR: could not determine sequence length of " << name
          << '\n'
          << "  FASTA file: " << d_->fastaFilename_;
        throw std::runtime_error{s.str()};
    }
    return len;
}
}  // namespace BAM
}  // namespace PacBio

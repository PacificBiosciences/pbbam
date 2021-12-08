#include "PbbamInternalConfig.h"

#include <pbbam/IndexedFastqReader.h>

#include <pbbam/BamRecord.h>
#include <pbbam/FaiIndex.h>
#include <pbbam/FormatUtils.h>
#include <pbbam/GenomicInterval.h>
#include "IndexedFastqBgzfReader.h"
#include "IndexedFastqReaderImpl.h"
#include "IndexedFastqTextReader.h"
#include "SequenceUtils.h"

#include <sstream>
#include <stdexcept>
#include <utility>

#include <cassert>

namespace PacBio {
namespace BAM {

namespace {

void ClipAndGapify(std::pair<std::string, Data::QualityValues>& seqQual, const Data::Cigar& cigar,
                   bool exciseSoftClips)
{
    auto& subseq = seqQual.first;
    auto subQual = seqQual.second.Fastq();

    const char nullQual = Data::QualityValue{0}.Fastq();

    size_t seqIndex = 0;
    for (const auto& op : cigar) {
        const auto type = op.Type();
        const auto opLength = op.Length();

        // do nothing for hard clips
        if (type == Data::CigarOperationType::HARD_CLIP) {
            continue;
        }

        // maybe remove soft clips
        if (type == Data::CigarOperationType::SOFT_CLIP) {
            if (!exciseSoftClips) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '-');
                subQual.reserve(subQual.size() + opLength);
                subQual.insert(seqIndex, opLength, nullQual);
                seqIndex += opLength;
            }
        }

        // for non-clipping operations
        else {
            // maybe add gaps/padding
            if (type == Data::CigarOperationType::INSERTION) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '-');
                subQual.reserve(subseq.size() + opLength);
                subQual.insert(seqIndex, opLength, nullQual);
            } else if (type == Data::CigarOperationType::PADDING) {
                subseq.reserve(subseq.size() + opLength);
                subseq.insert(seqIndex, opLength, '*');
                subQual.reserve(subseq.size() + opLength);
                subQual.insert(seqIndex, opLength, nullQual);
            }

            // update index
            seqIndex += opLength;
        }
    }

    seqQual.second = Data::QualityValues::FromFastq(subQual);
}

std::unique_ptr<IndexedFastqReaderImpl> MakeReaderImpl(std::string filename)
{
    // validate extension
    if (!FormatUtils::IsFastqFilename(filename)) {
        std::ostringstream s;
        s << "[pbbam] FASTQ reader ERROR: not a recognized FASTQ extension:\n"
          << "  file: " << filename;
        throw std::runtime_error{s.str()};
    }

    // determine subsequence "loader" from compression type: plain-text, bgzf, or unsupported
    const auto compressionType = FormatUtils::CompressionType(filename);
    switch (compressionType) {

        case HtslibCompression::NONE:
            return std::make_unique<IndexedFastqTextReader>(std::move(filename));
        case HtslibCompression::BGZIP:
            return std::make_unique<IndexedFastqBgzfReader>(std::move(filename));

        case HtslibCompression::GZIP: {
            std::ostringstream msg;
            msg << "[pbbam] FASTQ reader ERROR: random-access is not supported for plain gzipped "
                   "file "
                << filename << "\n\n"
                << "Compressed files must be bgzipped, with accompanying *.gzi "
                   "index.\n\n"
                << "To keep the original gzipped file unchanged:\n"
                << "  $ gunzip -c " << filename << " > <unzipped_file>\n"
                << "or discard the gzipped file:\n"
                << "  $ gunzip " << filename << '\n'
                << '\n'
                << "Re-compress & create *.gzi index:\n"
                << "  $ bgzip --index <unzipped_file>\n\n";
            throw std::runtime_error{msg.str()};
        }
        default: {
            assert(false);  // should never get here, the way htslib currently determines type
            std::ostringstream s;
            s << "[pbbam] FASTQ reader ERROR: could not determine compression type:\n"
              << "  file: " << filename;
            throw std::runtime_error{s.str()};
        }
    }
}

}  // namespace

IndexedFastqReader::IndexedFastqReader(std::string filename)
    : d_{MakeReaderImpl(std::move(filename))}
{
}

IndexedFastqReader::IndexedFastqReader(const IndexedFastqReader& other)
    : d_{MakeReaderImpl(other.d_->faiFilename_)}
{
}

IndexedFastqReader::IndexedFastqReader(IndexedFastqReader&&) noexcept = default;

IndexedFastqReader& IndexedFastqReader::operator=(const IndexedFastqReader& rhs)
{
    if (this != &rhs) {
        *this = IndexedFastqReader{rhs};
    }
    return *this;
}

IndexedFastqReader& IndexedFastqReader::operator=(IndexedFastqReader&&) noexcept = default;

IndexedFastqReader::~IndexedFastqReader() = default;

bool IndexedFastqReader::HasSequence(const std::string& name) const
{
    return d_->index_.HasEntry(name);
}

std::vector<std::string> IndexedFastqReader::Names() const { return d_->index_.Names(); }

std::string IndexedFastqReader::Name(const size_t idx) const { return d_->index_.Names().at(idx); }

int IndexedFastqReader::NumSequences() const { return d_->index_.Names().size(); }

std::pair<std::string, Data::QualityValues> IndexedFastqReader::ReferenceSubsequence(
    const BamRecord& bamRecord, const Data::Orientation orientation, const bool gapped,
    const bool exciseSoftClips)
{
    // fetch raw data for record's region
    auto seqQual = Subsequence(bamRecord.ReferenceName(), bamRecord.ReferenceStart(),
                               bamRecord.ReferenceEnd());

    // maybe clip/gapify
    if (bamRecord.Impl().IsMapped() && gapped) {
        ClipAndGapify(seqQual, bamRecord.Impl().CigarData(), exciseSoftClips);
    }

    // maybe reverse
    const auto reverse =
        (orientation != Data::Orientation::GENOMIC) && bamRecord.Impl().IsReverseStrand();
    if (reverse) {
        ReverseComplementCaseSens(seqQual.first);
        Reverse(seqQual.second);
    }

    return seqQual;
}

int IndexedFastqReader::SequenceLength(const std::string& name) const
{
    const auto& entry = d_->index_.Entry(name);
    return static_cast<int>(entry.Length);
}

std::pair<std::string, Data::QualityValues> IndexedFastqReader::Subsequence(const std::string& id,
                                                                            Data::Position start,
                                                                            Data::Position end)
{
    return d_->Subsequence(id, start, end);
}

std::pair<std::string, QualityValues> IndexedFastqReader::Subsequence(
    const Data::GenomicInterval& interval)
{
    return Subsequence(interval.Name(), interval.Start(), interval.Stop());
}

}  // namespace BAM
}  // namespace PacBio

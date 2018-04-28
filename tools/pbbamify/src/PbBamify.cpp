// Author: Ivan Sovic

#include "PbBamify.h"
#include <pbbam/../../src/SequenceUtils.h>
#include <pbbam/BamRecord.h>
#include <pbbam/Cigar.h>
#include <pbbam/MD5.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>
#include <ctime>
#include <istream>
#include <ostream>
#include <string>

namespace PacBio {
namespace BAM {
namespace pbbamify {

// Taken from BamRecord.cpp, since the implementation there
// is not public.
static inline bool ConsumesQuery(const CigarOperationType type)
{
    return (bam_cigar_type(static_cast<int>(type)) & 0x1) != 0;
}

// Taken from BamRecord.cpp, since the implementation there
// is not public.
static inline bool ConsumesReference(const CigarOperationType type)
{
    return (bam_cigar_type(static_cast<int>(type)) & 0x2) != 0;
}

PacBio::BAM::BamHeader Pbbamify::ComposeHeader(const PacBio::BAM::DataSet& dataset,
                                               PacBio::BAM::FastaReader& refReader,
                                               const PacBio::BAM::BamReader& input)
{

    PacBio::BAM::BamHeader retHeader;
    bool headerInitialized = false;

    // Merge all the read groups and additional PacBio info.
    const auto& bamFiles = dataset.BamFiles();
    for (auto& bamFile : bamFiles) {
        auto header = bamFile.Header();
        if (!headerInitialized) {
            retHeader = header.DeepCopy();
            headerInitialized = true;
        } else {
            retHeader += header;
        }
    }

    // Merge the alignment PG to the header.
    auto inputHeader = input.Header();
    for (auto& program : inputHeader.Programs()) {
        retHeader.AddProgram(program);
    }

    // Add the sequence info to the header.
    PacBio::BAM::FastaSequence record;
    while (refReader.GetNext(record)) {
        // Convert the sequence length to string,
        // as required by SequenceInfo.
        std::ostringstream ossLength;
        ossLength << record.Bases().size();

        // Clip on whitespace.
        std::istringstream issHeader(record.Name());
        std::string header;
        issHeader >> header;

        // Calculate the MD5 and append to retHeader.
        PacBio::BAM::SequenceInfo seq(header, ossLength.str());
        auto hash = PacBio::BAM::MD5Hash(record.Bases());
        seq.Checksum(hash);
        retHeader.AddSequence(seq);
    }

    return retHeader;
}

bool Pbbamify::IsHardClipped(const Cigar& cigarData)
{
    // If it's empty, just return.
    if (cigarData.size() == 0) {
        return false;
    }

    // If there is no hard clipping, just return.
    if (cigarData.front().Type() == CigarOperationType::HARD_CLIP ||
        cigarData.back().Type() == CigarOperationType::HARD_CLIP) {
        return true;
    }

    return false;
}

Cigar Pbbamify::ConvertHardToSoftClipping(const Cigar& cigarData)
{
    Cigar softCigar;

    // If it's empty, just return.
    if (cigarData.size() == 0) {
        return softCigar;
    }

    CigarOperationType prevOp = CigarOperationType::UNKNOWN_OP;

    for (const auto& cigar : cigarData) {
        // Change H to S.
        CigarOperationType op = (cigar.Type() == CigarOperationType::HARD_CLIP)
                                    ? CigarOperationType::SOFT_CLIP
                                    : cigar.Type();
        auto len = cigar.Length();

        // Merge or add.
        if (softCigar.size() > 0 && op == prevOp) {
            auto prevLen = softCigar.back().Length();
            softCigar.back() = CigarOperation(op, len + prevLen);
        } else {
            softCigar.emplace_back(CigarOperation(op, len));
        }

        prevOp = op;
    }

    return softCigar;
}

size_t Pbbamify::SequenceLengthFromCigar(const Cigar& cigarData)
{
    size_t len = 0;

    if (cigarData.size() == 0) {
        return len;
    }

    for (const auto& cigar : cigarData) {
        if (ConsumesQuery(cigar.Type()) || cigar.Type() == CigarOperationType::HARD_CLIP) {
            len += cigar.Length();
        }
    }

    return len;
}

bool Pbbamify::CheckIsCigarBasic(const Cigar& cigarData)
{
    for (const auto& cigar : cigarData) {
        if (cigar.Type() == CigarOperationType::ALIGNMENT_MATCH) {
            return true;
        }
    }
    return false;
}

/*
 * Takes the pre-calculated cigarData object so that it's
 * more efficient (it could always be obtained from the record
 * at any time).
*/
Cigar Pbbamify::BasicToExtendedCigar(const PacBio::BAM::IndexedFastaReader& indexedRefReader,
                                     const BamRecord& record, const Cigar& cigarData)
{
    Cigar extCigar;

    std::string qseq = record.Impl().Sequence();
    std::string rseq =
        indexedRefReader.ReferenceSubsequence(record, Orientation::GENOMIC, false, false);

    size_t qpos = 0, rpos = 0;  // The rpos should be 0 because the reference portion is yanked out.
    for (const auto& cigar : cigarData) {
        // This shouldn't happen, but let's keep it safe.
        if (cigar.Length() == 0) {
            continue;
        }

        if (cigar.Type() == CigarOperationType::ALIGNMENT_MATCH) {
            // Decode the prev op.
            CigarOperationType prevOp = (qseq[qpos] == rseq[rpos])
                                            ? CigarOperationType::SEQUENCE_MATCH
                                            : CigarOperationType::SEQUENCE_MISMATCH;
            size_t prevCount = 0;

            for (size_t i = 0; i < cigar.Length(); ++i) {
                // Decode the new op.
                CigarOperationType op = (qseq[qpos + i] == rseq[rpos + i])
                                            ? CigarOperationType::SEQUENCE_MATCH
                                            : CigarOperationType::SEQUENCE_MISMATCH;

                if (op == prevOp) {
                    ++prevCount;

                } else {
                    extCigar.emplace_back(CigarOperation(prevOp, prevCount));
                    prevOp = op;
                    prevCount = 1;
                }
            }

            // Add the last operation.
            extCigar.emplace_back(CigarOperation(prevOp, prevCount));

        } else {
            extCigar.emplace_back(cigar);
        }

        if (ConsumesQuery(cigar.Type())) {
            qpos += cigar.Length();
        }
        if (ConsumesReference(cigar.Type())) {
            rpos += cigar.Length();
        }
    }

    return extCigar;
}

bool Pbbamify::AugmentAlignments(const std::shared_ptr<QueryLookup> queryLookup,
                                 const PacBio::BAM::IndexedFastaReader& indexedRefReader,
                                 PacBio::BAM::BamReader& input, PacBio::BAM::BamWriter& writer,
                                 int32_t verboseLevel)
{

    // Clock is just for the verbose functionality.
    clock_t timerStart = clock();

    // Sets the frequency of the proof of life when
    // processing larger input BAMs.
    int32_t verboseFrequency =
        (verboseLevel <= 2)
            ? 1000000
            : (verboseLevel == 3)
                  ? 100000
                  : (verboseLevel == 4)
                        ? 10000
                        : (verboseLevel == 5)
                              ? 1000
                              : (verboseLevel == 6) ? 100 : (verboseLevel == 7) ? 10 : 1;

    // Counters for verbose output.
    size_t numRecords = 0, numWithoutSeq = 0;

    // Holder for the current record.
    BamRecord record;

    while (input.GetNext(record)) {
        ++numRecords;

        // Proof of life.
        if (verboseLevel > 1 && (numRecords % verboseFrequency) == 0) {
            double elapsedTime =
                static_cast<double>(clock() - timerStart) / (60.0 * CLOCKS_PER_SEC);
            elapsedTime = static_cast<int64_t>(elapsedTime * 100.0) / 100.0;
            std::cerr << "[INFO] Processed " << numRecords << " alignments in " << elapsedTime
                      << " min." << std::endl;
        }

        // Some mappers do not output sequences for secondary alignments.
        if (record.Impl().SequenceLength() == 0) {
            ++numWithoutSeq;
            continue;
        }

        // Update the BAM record with additional data from the PacBio dataset.
        int rv = AugmentAlignment(queryLookup, indexedRefReader, record, verboseLevel);

        // In case of failure, skip the alignment. Failures should be reported by AugmentAlignment.
        if (rv == false) {
            continue;
        }

        // Finally, write the output.
        writer.Write(record);
    }

    if (verboseLevel > 0 && numWithoutSeq) {
        std::cerr << "[Warning] Found " << numWithoutSeq
                  << " alignments without a seq field which were not converted (most likely "
                     "secondary alignments)."
                  << std::endl;
    }

    if (verboseLevel > 1) {
        double elapsedTime = static_cast<double>(clock() - timerStart) / (60.0 * CLOCKS_PER_SEC);
        elapsedTime = static_cast<int64_t>(elapsedTime * 100.0) / 100.0;
        std::cerr << "[INFO] Done processing " << numRecords << " alignments in " << elapsedTime
                  << " min." << std::endl;
    }

    return true;
}

bool Pbbamify::AugmentAlignment(const std::shared_ptr<QueryLookup> queryLookup,
                                const PacBio::BAM::IndexedFastaReader& indexedRefReader,
                                BamRecord& record, int32_t verboseLevel)
{

    // Find the BAM record in the original PacBio dataset.
    BamRecord datasetRecord;
    bool isFound = queryLookup->Find(record.FullName(), datasetRecord);
    if (isFound == 0) {
        if (verboseLevel > 0) {
            std::cerr << "[Warning] No records found for query '" << record.FullName()
                      << "'. Skipping." << std::endl;
        }
        return false;
    }

    // If it's not mapped, just output the original.
    if (record.IsMapped() == false) {
        record = datasetRecord;
        return true;
    }

    // Keep the cigar object since we'll reuse it. More efficient.
    auto cigar = record.Impl().CigarData();

    // Sanity check that the mapper did not produce something funky.
    size_t recordSeqLen = SequenceLengthFromCigar(cigar);
    if (recordSeqLen != datasetRecord.Impl().SequenceLength()) {
        if (verboseLevel > 0) {
            std::cerr << "[Warning] Sequence '" << record.FullName() << "' (length " << recordSeqLen
                      << ") is not of the same length as the PacBio BAM sequence (length "
                      << datasetRecord.Impl().SequenceLength() << ")! Skipping." << std::endl;
        }
        return false;
    }

    // Update the CIGAR only if necessary.
    if (CheckIsCigarBasic(cigar)) {
        cigar = BasicToExtendedCigar(indexedRefReader, record, cigar);
        record.Impl().CigarData(cigar);
    }

    // Stomp over any existing tags with matching IDs and add those
    // which do not yet exist in the aligned BAM. We consider the PacBio
    // dataset to be the correct answer to any of these. The rest are
    // produced by a mapper.
    // For example, BLASR will generate a RG tag even if the input was FASTA.
    for (auto& tag : datasetRecord.Impl().Tags()) {
        if (record.Impl().Tags().Contains(tag.first)) {
            record.Impl().EditTag(tag.first, tag.second);
        } else {
            record.Impl().AddTag(tag.first, tag.second);
        }
    }

    // Some downstream tools might not work well with the
    // "undefined" mapping quality value of 255. Here
    // we set it to a valid arbitrary value.
    if (record.Impl().MapQuality() == 255) {
        record.Impl().MapQuality(254);
    }

    // If the alignment has hard clipping, simply take both the seq and
    // qual fields from the dataset. This will stomp over any custom
    // qual values in the input BAM file.
    if (IsHardClipped(cigar)) {
        // Take the seq and qual fields from the dataset to override
        // any hard clippings induced by the mapper.
        std::string qseq = datasetRecord.Impl().Sequence();
        std::string quals = datasetRecord.Impl().Qualities().Fastq();

        // Reverse if needed.
        if (record.Impl().IsReverseStrand()) {
            PacBio::BAM::internal::ReverseComplement(qseq);
            std::reverse(quals.begin(), quals.end());
        }

        // PacBio datasets, when converted to SAM, contain '!' ASCII QVs.
        // In case QVs aren't provided otherwise, this block adds the '!' values.
        if (quals.size() == 0) {
            quals = std::string(qseq.size(), '!');
        }

        // Replace the seq and qual fields.
        record.Impl().SetSequenceAndQualities(qseq, quals);

        cigar = ConvertHardToSoftClipping(cigar);
        record.Impl().CigarData(cigar);

    } else {
        // PacBio datasets, when converted to SAM, contain '!' ASCII QVs.
        // In case QVs aren't provided otherwise, this block adds the '!' values.
        if (record.Impl().Qualities().size() == 0) {
            std::string qseq = record.Impl().Sequence();
            std::string quals = std::string(qseq.size(), '!');
            record.Impl().SetSequenceAndQualities(qseq, quals);
        }
    }

    return true;
}

}  // namespace pbbamify
}  // namespace BAM
}  // namespace PacBio

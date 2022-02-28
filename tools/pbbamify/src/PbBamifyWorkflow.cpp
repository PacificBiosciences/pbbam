#include "PbBamifyWorkflow.h"

#include <cstdint>
#include <ctime>

#include <sstream>
#include <string>

#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/DataSet.h>
#include <pbbam/FastaReader.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/MD5.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>

#include <pbcopper/data/Cigar.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/SequenceUtils.h>

#include "PbBamifySettings.h"
#include "PbBamifyVersion.h"

namespace PacBio {
namespace PbBamify {

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};

    // setup our @PG entry to add to header
    BAM::ProgramInfo pbbamifyProgram;
    pbbamifyProgram.Id(std::string{"pbbamify-"} + PbBamify::Version)
        .Name("pbbamify")
        .Version(PbBamify::Version);

    BAM::DataSet dataset{settings.PbbamFilename};
    BAM::BamReader inputBamReader{settings.InputFilename};
    BAM::BamHeader newHeader;

    {  // A separate block to close the reference file after the header is formed.
        // Using a sequential reader to construct the header SN lines in order, fast.
        BAM::FastaReader ref_reader{settings.ReferenceFilename};
        newHeader = ComposeHeader(dataset, ref_reader, inputBamReader);
    }

    auto queryLookup = std::make_shared<QueryLookup>(std::move(dataset));
    queryLookup->Load();

    {  // A block is used here to close the bamWriter and the reference reader.
        // (Even though this will be done as soon as the 'try' block ends, this safeguards if any
        // code should be added in between at some point.)
        BAM::IndexedFastaReader indexedRefReader{settings.ReferenceFilename};
        BAM::BamWriter bamWriter{settings.OutputFilename, newHeader};
        const bool result = AugmentAlignments(queryLookup, indexedRefReader, inputBamReader,
                                              bamWriter, settings.VerboseLevel);
        if (result == false) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

BAM::BamHeader Workflow::ComposeHeader(const BAM::DataSet& dataset, BAM::FastaReader& refReader,
                                       const BAM::BamReader& input)
{

    BAM::BamHeader retHeader;
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
    BAM::FastaSequence record;
    while (refReader.GetNext(record)) {

        // Convert the sequence length to string,
        // as required by SequenceInfo.
        std::ostringstream ossLength;
        ossLength << record.Bases().size();

        // Clip on whitespace.
        std::istringstream issHeader{record.Name()};
        std::string header;
        issHeader >> header;

        // Calculate the MD5 and append to retHeader.
        BAM::SequenceInfo seq{header, ossLength.str()};
        auto hash = BAM::MD5Hash(record.Bases());
        seq.Checksum(hash);
        retHeader.AddSequence(seq);
    }

    return retHeader;
}

bool Workflow::IsHardClipped(const Data::Cigar& cigarData)
{
    // If it's empty, just return.
    if (cigarData.size() == 0) {
        return false;
    }

    // If there is no hard clipping, just return.
    if (cigarData.front().Type() == Data::CigarOperationType::HARD_CLIP ||
        cigarData.back().Type() == Data::CigarOperationType::HARD_CLIP) {
        return true;
    }

    return false;
}

Data::Cigar Workflow::ConvertHardToSoftClipping(const Data::Cigar& cigarData)
{
    Data::Cigar softCigar;

    // If it's empty, just return.
    if (cigarData.size() == 0) {
        return softCigar;
    }

    Data::CigarOperationType prevOp = Data::CigarOperationType::UNKNOWN_OP;

    for (const auto& cigar : cigarData) {

        // Change H to S.
        Data::CigarOperationType op = (cigar.Type() == Data::CigarOperationType::HARD_CLIP)
                                          ? Data::CigarOperationType::SOFT_CLIP
                                          : cigar.Type();
        auto len = cigar.Length();

        // Merge or add.
        if (softCigar.size() > 0 && op == prevOp) {
            auto prevLen = softCigar.back().Length();
            softCigar.back() = Data::CigarOperation{op, len + prevLen};
        } else {
            softCigar.emplace_back(Data::CigarOperation{op, len});
        }

        prevOp = op;
    }

    return softCigar;
}

size_t Workflow::SequenceLengthFromCigar(const Data::Cigar& cigarData)
{
    size_t len = 0;

    if (cigarData.size() == 0) {
        return len;
    }

    for (const auto& cigar : cigarData) {
        if (Data::ConsumesQuery(cigar.Type()) ||
            cigar.Type() == Data::CigarOperationType::HARD_CLIP) {
            len += cigar.Length();
        }
    }

    return len;
}

bool Workflow::CheckIsCigarBasic(const Data::Cigar& cigarData)
{
    for (const auto& cigar : cigarData) {
        if (cigar.Type() == Data::CigarOperationType::ALIGNMENT_MATCH) {
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
Data::Cigar Workflow::BasicToExtendedCigar(const BAM::IndexedFastaReader& indexedRefReader,
                                           const BAM::BamRecord& record,
                                           const Data::Cigar& cigarData)
{
    Data::Cigar extCigar;

    std::string qseq{record.Impl().Sequence()};
    std::string rseq =
        indexedRefReader.ReferenceSubsequence(record, BAM::Orientation::GENOMIC, false, false);

    size_t qpos = 0, rpos = 0;  // The rpos should be 0 because the reference portion is yanked out.
    for (const auto& cigar : cigarData) {

        // This shouldn't happen, but let's keep it safe.
        if (cigar.Length() == 0) {
            continue;
        }

        if (cigar.Type() == Data::CigarOperationType::ALIGNMENT_MATCH) {
            // Decode the prev op.
            Data::CigarOperationType prevOp = (qseq[qpos] == rseq[rpos])
                                                  ? Data::CigarOperationType::SEQUENCE_MATCH
                                                  : Data::CigarOperationType::SEQUENCE_MISMATCH;
            uint32_t prevCount = 0;
            for (size_t i = 0; i < cigar.Length(); ++i) {

                // Decode the new op.
                Data::CigarOperationType op = (qseq[qpos + i] == rseq[rpos + i])
                                                  ? Data::CigarOperationType::SEQUENCE_MATCH
                                                  : Data::CigarOperationType::SEQUENCE_MISMATCH;

                if (op == prevOp) {
                    ++prevCount;
                } else {
                    extCigar.emplace_back(Data::CigarOperation{prevOp, prevCount});
                    prevOp = op;
                    prevCount = 1;
                }
            }

            // Add the last operation.
            extCigar.emplace_back(Data::CigarOperation{prevOp, prevCount});
        } else {
            extCigar.emplace_back(cigar);
        }

        if (Data::ConsumesQuery(cigar.Type())) {
            qpos += cigar.Length();
        }
        if (Data::ConsumesReference(cigar.Type())) {
            rpos += cigar.Length();
        }
    }

    return extCigar;
}

bool Workflow::AugmentAlignments(const std::shared_ptr<QueryLookup> queryLookup,
                                 const BAM::IndexedFastaReader& indexedRefReader,
                                 BAM::BamReader& input, BAM::BamWriter& writer,
                                 int32_t verboseLevel)
{

    // Clock is just for the verbose functionality.
    clock_t timerStart = clock();

    // Sets the frequency of the proof of life when
    // processing larger input BAMs.
    int32_t verboseFrequency = (verboseLevel <= 2)   ? 1000000
                               : (verboseLevel == 3) ? 100000
                               : (verboseLevel == 4) ? 10000
                               : (verboseLevel == 5) ? 1000
                               : (verboseLevel == 6) ? 100
                               : (verboseLevel == 7) ? 10
                                                     : 1;

    // Counters for verbose output.
    size_t numRecords = 0, numWithoutSeq = 0;

    // Holder for the current record.
    BAM::BamRecord record;
    while (input.GetNext(record)) {
        ++numRecords;

        // Proof of life.
        if (verboseLevel > 1 && (numRecords % verboseFrequency) == 0) {
            double elapsedTime =
                static_cast<double>(clock() - timerStart) / (60.0 * CLOCKS_PER_SEC);
            elapsedTime = static_cast<int64_t>(elapsedTime * 100.0) / 100.0;
            PBLOG_INFO << "Processed " << numRecords << " alignments in " << elapsedTime << " min.";
        }

        // Some mappers do not output sequences for secondary alignments.
        if (record.Impl().SequenceLength() == 0) {
            ++numWithoutSeq;
            continue;
        }

        // Update the BAM record with additional data from the PacBio dataset.
        // In case of failure, skip the alignment. Failures should be reported by AugmentAlignment.
        const bool rv = AugmentAlignment(queryLookup, indexedRefReader, record, verboseLevel);
        if (rv == false) {
            continue;
        }

        // Finally, write the output.
        writer.Write(record);
    }

    if (verboseLevel > 0 && numWithoutSeq) {
        PBLOG_WARN << "Found " << numWithoutSeq
                   << " alignments without a seq field which were not converted (most likely "
                      "secondary alignments).";
    }

    if (verboseLevel > 1) {
        double elapsedTime = static_cast<double>(clock() - timerStart) / (60.0 * CLOCKS_PER_SEC);
        elapsedTime = static_cast<int64_t>(elapsedTime * 100.0) / 100.0;
        PBLOG_INFO << "Done processing " << numRecords << " alignments in " << elapsedTime
                   << " min.";
    }

    return true;
}

bool Workflow::AugmentAlignment(const std::shared_ptr<QueryLookup> queryLookup,
                                const BAM::IndexedFastaReader& indexedRefReader,
                                BAM::BamRecord& record, int32_t verboseLevel)
{

    // Find the BAM record in the original PacBio dataset.
    BAM::BamRecord datasetRecord;
    const bool isFound = queryLookup->Find(record.FullName(), datasetRecord);
    if (!isFound) {
        if (verboseLevel > 0) {
            PBLOG_WARN << "No records found for query '" << record.FullName() << "'. Skipping.";
        }
        return false;
    }

    // If it's not mapped, just output the original.
    if (!record.IsMapped()) {
        record = datasetRecord;
        return true;
    }

    // Keep the cigar object since we'll reuse it. More efficient.
    auto cigar = record.Impl().CigarData();

    // Sanity check that the mapper did not produce something funky.
    const size_t recordSeqLen = SequenceLengthFromCigar(cigar);
    if (recordSeqLen != datasetRecord.Impl().SequenceLength()) {
        if (verboseLevel > 0) {
            PBLOG_WARN << "Sequence '" << record.FullName() << "' (length " << recordSeqLen
                       << ") is not of the same length as the PacBio BAM sequence (length "
                       << datasetRecord.Impl().SequenceLength() << ")! Skipping.";
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
        std::string qseq{datasetRecord.Impl().Sequence()};
        std::string quals{datasetRecord.Impl().Qualities().Fastq()};

        // Reverse if needed.
        if (record.Impl().IsReverseStrand()) {
            Utility::ReverseComplement(qseq);
            std::reverse(quals.begin(), quals.end());
        }

        // PacBio datasets, when converted to SAM, contain '!' ASCII QVs.
        // In case QVs aren't provided otherwise, this block adds the '!' values.
        if (quals.size() == 0) {
            quals = std::string(qseq.size(), '!');
        }

        // Replace the seq, qual, & cigar fields.
        record.Impl().SetSequenceAndQualities(qseq, quals);
        cigar = ConvertHardToSoftClipping(cigar);
        record.Impl().CigarData(cigar);

    } else {
        // PacBio datasets, when converted to SAM, contain '!' ASCII QVs.
        // In case QVs aren't provided otherwise, this block adds the '!' values.
        if (record.Impl().Qualities().size() == 0) {
            std::string qseq{record.Impl().Sequence()};
            std::string quals(qseq.size(), '!');
            record.Impl().SetSequenceAndQualities(qseq, quals);
        }
    }

    return true;
}

}  // namespace PbBamify
}  // namespace PacBio

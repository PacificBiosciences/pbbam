#ifndef PBBAMIFY_WORKFLOW_H
#define PBBAMIFY_WORKFLOW_H

#include <pbcopper/cli2/Results.h>

#include "QueryLookup.h"

namespace PacBio {
namespace BAM {

class BamHeader;
class BamReader;
class BamRecord;
class BamWriter;
class DataSet;
class FastaReader;
class IndexedFastaReader;

}  // namespace BAM

namespace Data {

class Cigar;

}  // namespace Data

namespace PbBamify {

struct Workflow
{
    ///
    /// \brief Takes a PacBio dataset, a reference file and an input arbitrary
    ///        aligned BAM. Produces a new PacBio-compatible aligned BAM.
    ///
    /// \throws std::runtime_error
    ///
    static int Runner(const CLI_v2::Results& args);

    ///
    /// \brief Merges all the headers from the dataset and the input, adds the
    ///        SQ fields with lengths and MD5 checksums.
    ///
    /// \returns A BAM header which is composed of: merged headers from BAMs in
    ///          the dataset, ProgramInfo from the input BAM, and SQ lines
    ///          formed from the refReader object (together with their length and
    ///          MD5 checksum).
    ///
    static BAM::BamHeader ComposeHeader(const BAM::DataSet& dataset, BAM::FastaReader& refReader,
                                        const BAM::BamReader& input);

    ///
    /// \brief Converts a set of generic BAM records into a PacBio compatible
    ///        BAM by calling AugmentAlignment for each BAM record in the input
    ///        BAM file. If a BAM record was not mapped, then the original record
    ///        from the dataset will be set to `record`.
    ///
    /// \returns true if the record was successfully augmented, false otherwise.
    ///
    static bool AugmentAlignments(const std::shared_ptr<QueryLookup> queryLookup,
                                  const BAM::IndexedFastaReader& indexedRefReader,
                                  BAM::BamReader& input, BAM::BamWriter& writer,
                                  int32_t verboseLevel);

    ///
    /// \brief Converts a generic BAM record into a PacBio compatible BAM by:
    ///        adding tags from the PacBio dataset, replacing the read group,
    ///        clipping the tags if needed, converting the CIGAR from basic to
    ///        extended format if needed, changing the mapq from 255 to another
    ///        value to avoid potential downstream issues, etc.
    ///
    /// \returns true if the record was successfully augmented, false otherwise.
    ///
    static bool AugmentAlignment(const std::shared_ptr<QueryLookup> queryLookup,
                                 const BAM::IndexedFastaReader& indexedRefReader,
                                 BAM::BamRecord& record, int32_t verboseLevel);

    ///
    /// \brief Checks whether the alignment was hard clipped.
    ///
    /// \returns true if the front or back CIGAR op is 'H', false otherwise.
    ///
    static bool IsHardClipped(const Data::Cigar& cigarData);

    ///
    /// \brief If the CIGAR string contains hard clipping operation at the beginning
    ///        or end of the cigarData vector, these are turned to soft clips and
    ///        merged with any potential existin soft clipping operations.
    ///
    /// \returns a new CIGAR string with only soft clipped bases.
    ///
    static Data::Cigar ConvertHardToSoftClipping(const Data::Cigar& cigarData);

    ///
    /// \brief Calculates the total sequence length from CIGAR (including clipping),
    ///        and not just the aligned length. This is used for sanity checking
    ///        the input BAM records.
    ///
    /// \returns The length of the query sequence calculated from the CIGAR string.
    ///
    static size_t SequenceLengthFromCigar(const Data::Cigar& cigarData);

    ///
    /// \brief Linear pass over the Cigar operations to see if there are any 'M' ops.
    ///
    /// \returns true if there are 'M' operations in the CIGAR object.
    ///
    static bool CheckIsCigarBasic(const Data::Cigar& cigarData);

    ///
    /// \brief Takes the index and a BAM record, and creates a new Cigar object
    ///        with extended CIGAR operations ('=' and 'X' instead of 'M').
    ///
    /// \returns A new Cigar object with '=' and 'X' operations instead of 'M's.
    ///
    static Data::Cigar BasicToExtendedCigar(const BAM::IndexedFastaReader& indexedRefReader,
                                            const BAM::BamRecord& record,
                                            const Data::Cigar& cigarData);
};

}  // namespace PbBamify
}  // namespace PacBio

#endif  // PBBAMIFY_WORKFLOW_H

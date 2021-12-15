#ifndef PBBAM_INDEXEDFASTAREADER_H
#define PBBAM_INDEXEDFASTAREADER_H

#include <pbbam/Config.h>

#include <pbbam/BamRecord.h>
#include <pbbam/Orientation.h>
#include <pbbam/Position.h>

#include <pbcopper/data/GenomicInterval.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace BAM {

/// \brief The IndexedFastaReader class provides random-access to FASTA file
///        data.
///
class IndexedFastaReader
{

public:
    /// \name Constructors & Related Methods
    /// \{

    explicit IndexedFastaReader(std::string filename);

    IndexedFastaReader(const IndexedFastaReader&);
    IndexedFastaReader(IndexedFastaReader&&) noexcept;
    IndexedFastaReader& operator=(const IndexedFastaReader&);
    IndexedFastaReader& operator=(IndexedFastaReader&&) noexcept;
    ~IndexedFastaReader();

    /// \}

public:
    /// name Sequence Access
    /// \{

    /// \brief Fetches FASTA sequence for desired interval.
    ///
    /// \param[in] id       reference sequence name
    /// \param[in] begin    start position
    /// \param[in] end      end position
    ///
    /// \returns sequence string at desired interval
    ///
    /// \throws std::runtime_error on failure to fetch sequence
    ///
    std::string Subsequence(const std::string& id, Data::Position begin, Data::Position end) const;

    /// \brief Fetches FASTA sequence for desired interval.
    ///
    /// \param[in] interval desired interval
    ///
    /// \returns sequence string at desired interval
    ///
    /// \throws std::runtime_error on failure to fetch sequence
    ///
    std::string Subsequence(const Data::GenomicInterval& interval) const;

    /// \brief Fetches FASTA sequence for desired interval.
    ///
    /// \param[in] htslibRegion htslib/samtools-formatted REGION string
    ///                         representing the desired interval
    ///
    /// \returns sequence string at desired interval
    ///
    /// \throws std::runtime_error on failure to fetch sequence
    ///
    std::string Subsequence(const char* htslibRegion) const;

    /// \brief Fetches FASTA sequence corresponding to a BamRecord, oriented and
    ///        gapped as requested.
    ///
    /// For example, "native" orientation and "gapped" will return the reference
    /// sequence with gaps inserted, as would align against the read in "native"
    /// orientation.
    ///
    /// \param[in] bamRecord        input BamRecord to derive interval/CIGAR
    ///                             data
    /// \param[in] orientation      orientation of output
    /// \param[in] gapped           if true, gaps/padding will be inserted, per
    ///                             record's CIGAR info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    /// \returns sequence string over the record's interval
    ///
    /// \throws std::runtime_error on failure to fetch sequence
    ///
    std::string ReferenceSubsequence(const BamRecord& bamRecord,
                                     Data::Orientation orientation = Data::Orientation::GENOMIC,
                                     bool gapped = false, bool exciseSoftClips = false) const;

    /// \}

public:
    /// \name File Attributes
    /// \{

    /// \returns true if FASTA file contains a sequence matching \p name
    bool HasSequence(const std::string& name) const;

    /// \returns the names of the sequence at a specific index in the FASTA file
    std::string Name(size_t idx) const;

    /// \returns the names of all sequences stored in the FASTA file
    std::vector<std::string> Names() const;

    /// \returns number of sequences stored in FASTA file
    int NumSequences() const;

    /// \returns length of FASTA sequence
    ///
    /// \throws std::runtime_error if length could not be determined
    ///
    int SequenceLength(const std::string& name) const;

    /// \}

private:
    class IndexedFastaReaderPrivate;
    std::unique_ptr<IndexedFastaReaderPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_INDEXEDFASTAREADER_H

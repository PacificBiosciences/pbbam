// File Description
/// \file IndexedFastaReader.h
/// \brief Defines the IndexedFastaReader class.
//
// Author: David Alexander

#ifndef INDEXEDFASTAREADER_H
#define INDEXEDFASTAREADER_H

#include <htslib/faidx.h>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "pbbam/Orientation.h"
#include "pbbam/Position.h"

namespace PacBio {
namespace BAM {

class GenomicInterval;
class BamRecord;

/// \brief The IndexedFastaReader class provides random-access to FASTA file
///        data.
///
class IndexedFastaReader
{

public:
    /// \name Constructors & Related Methods
    /// \{

    IndexedFastaReader() = delete;
    IndexedFastaReader(const std::string& filename);
    IndexedFastaReader(const IndexedFastaReader& src);
    IndexedFastaReader(IndexedFastaReader&&) = default;
    IndexedFastaReader& operator=(const IndexedFastaReader& rhs);
    IndexedFastaReader& operator=(IndexedFastaReader&&) = default;
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
    std::string Subsequence(const std::string& id, Position begin, Position end) const;

    /// \brief Fetches FASTA sequence for desired interval.
    ///
    /// \param[in] interval desired interval
    ///
    /// \returns sequence string at desired interval
    ///
    /// \throws std::runtime_error on failure to fetch sequence
    ///
    std::string Subsequence(const GenomicInterval& interval) const;

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
                                     const Orientation orientation = Orientation::GENOMIC,
                                     const bool gapped = false,
                                     const bool exciseSoftClips = false) const;

    /// \}

public:
    /// \name File Attributes
    /// \{

    /// \returns true if FASTA file contains a sequence matching \p name
    bool HasSequence(const std::string& name) const;

    /// \returns the names of the sequence at a specific index in the FASTA file
    std::string Name(const size_t idx) const;

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
    std::string filename_;
    faidx_t* handle_;

private:
    void Close();
    bool Open(std::string filename);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INDEXEDFASTAREADER_H

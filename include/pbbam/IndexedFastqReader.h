// File Description
/// \file IndexedFastqReader.h
/// \brief Defines the IndexedFastqReader class.
//
// Author: Derek Barnett

#ifndef INDEXEDFASTQREADER_H
#define INDEXEDFASTQREADER_H

#include "pbbam/Config.h"

#include <cstddef>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "pbbam/FastqReader.h"
#include "pbbam/Orientation.h"

#include "internal/QueryBase.h"

namespace PacBio {
namespace BAM {

class BamRecord;
class GenomicInterval;
class IndexedFastqReaderImpl;

/// \brief The IndexedFastaReader class provides random-access to FASTQ file
///        data.
///
class IndexedFastqReader
{

public:
    /// \name Constructors & Related Methods
    /// \{

    explicit IndexedFastqReader(std::string filename);

    IndexedFastqReader(const IndexedFastqReader&);
    IndexedFastqReader(IndexedFastqReader&&) noexcept;
    IndexedFastqReader& operator=(const IndexedFastqReader& rhs);
    IndexedFastqReader& operator=(IndexedFastqReader&&) noexcept;
    ~IndexedFastqReader();

    /// \}

public:
    /// name Sequence Access
    /// \{

    /// \brief Fetches sequence & qualities for desired interval.
    ///
    /// \param[in] id       reference sequence name
    /// \param[in] start    start position
    /// \param[in] end      end position
    ///
    /// \returns sequence/QV pair for desired interval
    ///
    /// \throws std::runtime_error on failure to fetch data
    ///
    std::pair<std::string, Data::QualityValues> Subsequence(const std::string& id,
                                                            Data::Position start,
                                                            Data::Position end);

    /// \brief Fetches sequence & qualities for desired interval.
    ///
    /// \param[in] interval desired interval
    ///
    /// \returns sequence/QV pair for desired interval
    ///
    /// \throws std::runtime_error on failure to fetch data
    ///
    std::pair<std::string, Data::QualityValues> Subsequence(const GenomicInterval& interval);

    /// \brief Fetches sequence & qualities sequence corresponding to a BamRecord, oriented and
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
    /// \returns sequence/QV pair over the record's interval
    ///
    /// \throws std::runtime_error on failure to fetch data
    ///
    std::pair<std::string, Data::QualityValues> ReferenceSubsequence(
        const BamRecord& bamRecord, const Orientation orientation = Orientation::GENOMIC,
        const bool gapped = false, const bool exciseSoftClips = false);

    /// \}

public:
    /// \name File Attributes
    /// \{

    /// \returns true if FASTQ file contains a sequence matching \p name
    bool HasSequence(const std::string& name) const;

    /// \returns the names of the sequence at a specific index in the FASTQ file
    std::string Name(const size_t idx) const;

    /// \returns the names of all sequences stored in the FASTQ file
    std::vector<std::string> Names() const;

    /// \returns number of sequences stored in FASTQ file
    int NumSequences() const;

    /// \returns length of FASTQ sequence
    ///
    /// \throws std::runtime_error if length could not be determined
    ///
    int SequenceLength(const std::string& name) const;

    /// \}

private:
    std::unique_ptr<IndexedFastqReaderImpl> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INDEXEDFASTQREADER_H

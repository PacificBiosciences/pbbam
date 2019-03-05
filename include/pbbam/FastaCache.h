// File Description
/// \file FastaCache.h
/// \brief Defines the FastaCache
//
// Author: Derek Barnett

#ifndef FASTACACHE_H
#define FASTACACHE_H

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "pbbam/Config.h"
#include "pbbam/FastaSequence.h"

namespace PacBio {
namespace BAM {

class FastaCacheData
{
public:
    explicit FastaCacheData(const std::string& filename);

    /// \brief Fetches FASTA sequence for desired interval.
    ///
    /// \param[in] name     reference sequence name
    /// \param[in] begin    start position
    /// \param[in] end      end position
    ///
    /// \returns sequence string at desired interval
    ///
    /// \throws std::runtime_error on failure to fetch sequence
    ///
    std::string Subsequence(const std::string& name, size_t begin, size_t end) const;

    /// \returns the names of all sequences stored in the FASTA file
    std::vector<std::string> Names() const;

    /// \returns length of FASTA sequence
    ///
    /// \throws std::runtime_error if name is unknown
    ///
    size_t SequenceLength(const std::string& name) const;

private:
    std::vector<FastaSequence> cache_;
    std::unordered_map<std::string, size_t> lookup_;
};

using FastaCache = std::shared_ptr<FastaCacheData>;

FastaCache MakeFastaCache(const std::string& filename);

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTACACHE_H
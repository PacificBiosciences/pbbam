#ifndef PBBAM_FASTACACHE_H
#define PBBAM_FASTACACHE_H

#include <pbbam/Config.h>

#include <pbbam/FastaSequence.h>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace PacBio {
namespace BAM {

class FastaCacheData
{
public:
    explicit FastaCacheData(const std::string& filename);

    ///
    /// \brief Check FASTA data for valid bases.
    ///
    /// Check to see if sequences contain only [ACGTNacgtn] characters.
    ///
    /// \returns a pair consisting of a bool denoting success/fail and, on failure,
    ///          the name of the first failing FASTA entry.
    ///
    std::pair<bool, std::string> Check() const;

    ///
    /// \brief Check FASTA data using custom callback.
    ///
    /// Uses the callable provided for validation, which should return
    /// success or fail for a given FASTA entry.
    ///
    /// \returns a pair consisting of a bool denoting success/fail and, on failure,
    ///          the name of the first failing FASTA entry.
    ///
    std::pair<bool, std::string> Check(
        const std::function<bool(const FastaSequence&)>& callback) const;

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

#endif  // PBBAM_FASTACACHE_H

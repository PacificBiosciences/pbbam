#include "PbbamInternalConfig.h"

#include "pbbam/FastaCache.h"

#include <sstream>
#include <stdexcept>

#include "pbbam/FastaReader.h"

namespace PacBio {
namespace BAM {

FastaCacheData::FastaCacheData(const std::string& filename) : cache_{FastaReader::ReadAll(filename)}
{
    for (size_t i = 0; i < cache_.size(); ++i)
        lookup_.emplace(cache_[i].Name(), i);
}

std::pair<bool, std::string> FastaCacheData::Check() const
{
    return Check([](const FastaSequence& seq) {
        const auto& bases = seq.Bases();
        const auto invalid = bases.find_first_not_of("ACGTNacgtn\n");
        return invalid == std::string::npos;
    });
}

std::pair<bool, std::string> FastaCacheData::Check(
    const std::function<bool(const FastaSequence&)>& callback) const
{
    for (const auto& seq : cache_) {
        if (!callback(seq)) {
            return {false, seq.Name()};
        }
    }
    return {true, ""};
}

std::string FastaCacheData::Subsequence(const std::string& name, size_t begin, size_t end) const
{
    const auto found = lookup_.find(name);
    if (found == lookup_.cend()) {
        std::ostringstream s;
        s << "[pbbam] FASTA sequence cache ERROR: could not retrieve subsequence, reference '"
          << name << "' not found";
        throw std::runtime_error{s.str()};
    }
    const std::string& seq = cache_[found->second].Bases();

    if (begin > end) {
        std::ostringstream s;
        s << "[pbbam] FASTA sequence cache ERROR: could not retrieve subsequence, requested begin "
             "> end\n"
          << "  begin: " << begin << '\n'
          << "  end: " << end << '\n'
          << "  rname: " << name;
        throw std::runtime_error{s.str()};
    }
    const size_t length = end - begin;
    return seq.substr(begin, length);
}

std::vector<std::string> FastaCacheData::Names() const
{
    std::vector<std::string> result;
    result.reserve(cache_.size());
    for (const auto& seq : cache_)
        result.push_back(seq.Name());
    return result;
}

size_t FastaCacheData::SequenceLength(const std::string& name) const
{
    const auto found = lookup_.find(name);
    if (found == lookup_.cend()) {
        std::ostringstream s;
        s << "[pbbam] FASTA sequence cache ERROR: could not retrieve sequence loength, reference '"
          << name << "' not found";
        throw std::runtime_error{s.str()};
    }
    return cache_[found->second].Bases().size();
}

FastaCache MakeFastaCache(const std::string& filename)
{
    return std::make_shared<FastaCacheData>(filename);
}

}  // namespace BAM
}  // namespace PacBio

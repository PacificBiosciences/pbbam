#include "pbbam/FastaCache.h"

#include <stdexcept>

#include "pbbam/FastaReader.h"

namespace PacBio {
namespace BAM {

FastaCacheData::FastaCacheData(const std::string& filename) : cache_{FastaReader::ReadAll(filename)}
{
    for (size_t i = 0; i < cache_.size(); ++i)
        lookup_.emplace(cache_[i].Name(), i);
}

std::string FastaCacheData::Subsequence(const std::string& name, size_t begin, size_t end) const
{
    const auto found = lookup_.find(name);
    if (found == lookup_.cend()) {
        throw std::runtime_error{""};
    }
    const std::string& seq = cache_[found->second].Bases();

    if (begin > end) throw std::runtime_error{""};
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
        throw std::runtime_error{""};
    }
    return cache_[found->second].Bases().size();
}

FastaCache MakeFastaCache(const std::string& filename)
{
    return std::make_shared<FastaCacheData>(filename);
}

}  // namespace BAM
}  // namespace PacBio

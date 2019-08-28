#include "PbbamInternalConfig.h"

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
        std::string msg = "Could not find '";
        msg += name;
        msg += "' in FastaCacheData::Subsequence()";
        throw std::runtime_error{msg};
    }
    const std::string& seq = cache_[found->second].Bases();

    if (begin > end) throw std::runtime_error{"begin > end in FastaCacheData::Subsequence"};
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
        std::string msg = "Could not find '";
        msg += name;
        msg += "' in FastaCacheData::SequenceLength()";
        throw std::runtime_error{msg};
    }
    return cache_[found->second].Bases().size();
}

FastaCache MakeFastaCache(const std::string& filename)
{
    return std::make_shared<FastaCacheData>(filename);
}

}  // namespace BAM
}  // namespace PacBio

// Author: Ivan Sovic

#ifndef PBBAMIFY_QUERYLOOKUP_H
#define PBBAMIFY_QUERYLOOKUP_H

#include <cstdint>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>

namespace PacBio {
namespace PbBamify {

///
/// \brief A simple container to hold the location of a read.
///
struct QueryLocation
{
    uint16_t fileNumber = 0;
    int64_t fileOffset = 0;
};

///
/// \brief QueryLookup parses all reads from PacBio indexes and creates a
///        hash lookup where the key is the read's qname, and the value is a
///        QueryLocation object pointing to the exact location of the read. The BAM
///        record can then be loaded by setting the virtual offset and calling GetNext().
///
class QueryLookup
{
public:
    explicit QueryLookup(BAM::DataSet dataset);
    QueryLookup(const QueryLookup&) = delete;
    QueryLookup& operator=(const QueryLookup&) = delete;

    ///
    /// \brief  Load() performs the work of setting up the BamReaders and constructing
    ///         the hash table lookup.
    ///
    /// \throws std::runtime_error if there are more than 1 record for a given qname.
    ///
    void Load();

    ///
    /// \brief Find(...) attempts to find a given qName in the lookup and return
    ///        the related BAM record. If it cannot be found, the function returns false.
    ///
    /// \returns true if the record was found and loaded, false otherwise.
    ///
    bool Find(const std::string& qName, BAM::BamRecord& record) const;

private:
    BAM::DataSet dataset_;
    std::vector<std::shared_ptr<BAM::BamReader>> readers_;
    std::unordered_map<std::string, QueryLocation> lookup_;
};

}  // namespace PbBamify
}  // namespace PacBio

#endif  // PBBAMIFY_QUERYLOOKUP_H

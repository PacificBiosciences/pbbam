// Copyright (c) 2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Ivan Sovic

#ifndef PBBAMIFY_SRC_QUERY_LOOKUP_H_
#define PBBAMIFY_SRC_QUERY_LOOKUP_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>

namespace PacBio {
namespace BAM {
namespace pbbamify {

class QueryLocation;
class QueryLookup;

/// \brief A factory function for the QueryLookup objects.
///
/// \returns A new QueryLookup object wrapped in a unique_ptr.
///
std::unique_ptr<QueryLookup> CreateQueryLookup(const PacBio::BAM::DataSet& dataset);

/// \brief A simple container to hold the location of a read.
class QueryLocation {
public:
    QueryLocation()
                 : fileNumber{0}, fileOffset{0} { }
    QueryLocation(uint16_t _fileNumber, int64_t _fileOffset)
                 : fileNumber{_fileNumber}, fileOffset{_fileOffset} { }

    uint16_t fileNumber;
    int64_t fileOffset;
};

/// \brief QueryLookup parses all reads from PacBio indexes and creates a
///        hash lookup where the key is the read's qname, and the value is a
///        QueryLocation object pointing to the exact location of the read. The BAM
///        record can then be loaded by setting the virtual offset and calling GetNext().
class QueryLookup {
public:
    friend std::unique_ptr<QueryLookup> CreateQueryLookup(const PacBio::BAM::DataSet& dataset);

    ~QueryLookup() = default;

    /// \brief  Load() performs the work of setting up the BamReaders and constructing
    ///         the hash table lookup.
    ///
    /// \throws std::runtime_error if there are more than 1 record for a given qname.
    void Load();

    /// \brief Find(...) attempts to find a given qName in the lookup and return
    ///        the related BAM record. If it cannot be found, the function returns false.
    ///
    /// \returns true if the record was found and loaded, false otherwise.
    bool Find(const std::string& qName, BamRecord& record) const;

private:
    QueryLookup(const QueryLookup&) = delete;
    QueryLookup& operator=(const QueryLookup&) = delete;

    /// \brief The constructor simply initializes a private reference to the dataset. No work is performed here.
    QueryLookup(const PacBio::BAM::DataSet& dataset);

    const PacBio::BAM::DataSet& dataset_;
    std::vector<std::shared_ptr<PacBio::BAM::BamReader>> readers_;
    std::unordered_map<std::string, QueryLocation> lookup_;
};

} // namespace pbbamify
} // namespace BAM
} // namespace PacBio

#endif

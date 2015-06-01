// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

// Author: Derek Barnett

#ifndef SUBDATASETS_H
#define SUBDATASETS_H

#include "pbbam/Config.h"
#include "pbbam/dataset/Filters.h"
#include "pbbam/internal/DataSetElement.h"
#include "pbbam/internal/DataSetListElement.h"
#include <string>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT SubDataSet : public internal::DataSetElement
{
public:
    SubDataSet(void);
    using DataSetElement::DataSetElement;

public:
    std::string CreatedAt(void) const;
    std::string Name(void) const;
    std::string Tags(void) const;
    std::string UniqueId(void) const;
    std::string Version(void) const;

    const Filters& FilterList(void) const;
    Filters& FilterList(void);

    void AddFilter(const Filter& filter);
    void RemoveFilter(const Filter& filter);

public:
    SubDataSet& CreatedAt(const std::string& timestamp);
    SubDataSet& Name(const std::string& name);
    SubDataSet& Tags(const std::string& tags);
    SubDataSet& UniqueId(const std::string& uuid);
    SubDataSet& Version(const std::string& version);
};

class PBBAM_EXPORT SubDataSets : public internal::DataSetListElement<SubDataSet>
{
public:
    SubDataSets(void);
    using DataSetListElement::DataSetListElement;

public:
    void AddSubDataSet(const SubDataSet& subdataset);
    void RemoveSubDataSet(const SubDataSet& subdataset);
};

} // namespace BAM
} // namespace PacBio

#endif // SUBDATASETS_H

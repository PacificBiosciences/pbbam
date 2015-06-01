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

#ifndef DATASET_H
#define DATASET_H

#include "pbbam/BamFile.h"
#include "pbbam/dataset/AlignmentSet.h"
#include "pbbam/dataset/BarcodeSet.h"
#include "pbbam/dataset/CcsReadSet.h"
#include "pbbam/dataset/ContigSet.h"
#include "pbbam/dataset/ReferenceSet.h"
#include "pbbam/dataset/SubreadSet.h"

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT DataSet : public DataSetBase
{

public:
    DataSet(void);
    DataSet(const BamFile& bamFile);
    using DataSetBase::DataSetBase;

public:
    /// \name Generic Dataset Metadata
    /// \{

    /// \returns dataset's metadata
    ///
    const DataSetMetadata& Metadata(void) const;

public:
    /// \returns editable metadata object
    DataSetMetadata& Metadata(void);

    /// \}

public:
    /// \name DataSet Conversion
    /// \{

    AlignmentSet ToAlignmentSet(void) const;
    BarcodeSet   ToBarcodeSet(void) const;
    CcsReadSet   ToCcsReadSet(void) const;
    ContigSet    ToContigSet(void) const;
    ReferenceSet ToReferenceSet(void) const;
    SubreadSet   ToSubreadSet(void) const;

    /// \}

public:
    /// \n BAM-Specific Methods
    /// \{

    /// \returns list of BamFile objects for BAM files in external data references
    ///
    std::vector<BamFile> BamFiles(void) const;

    /// \}

public:
    /// \name Merging Datasets
    /// \{

    /// Merges \p other dataset with this dataset.
    ///
    /// Checks for same metadata (for now).
    ///
    /// \returns reference to this dataset
    ///
    DataSet& operator+=(const DataSet& other);

    /// \}
};

} // namespace BAM
} // namespace PacBio

#endif // DATASET_H

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
//
// Author: Derek Barnett

#ifndef PBIINDEX_H
#define PBIINDEX_H

#include "pbbam/Config.h"
#include "pbbam/PbiFile.h"
#include "pbbam/PbiLookupData.h"
#include <memory>
#include <string>

namespace PacBio {
namespace BAM {

namespace internal { class PbiIndexPrivate; }

class PBBAM_EXPORT PbiIndex
{
public:
    /// \name Constructors & Related Methods
    /// \{

    /// Constructs a PbiIndex lookup structure from a PBI file.
    PbiIndex(const std::string& pbiFilename);

    PbiIndex(const PbiIndex& other);
    PbiIndex(PbiIndex&& other);
    PbiIndex& operator=(const PbiIndex& other);
    PbiIndex& operator=(PbiIndex&& other);
    ~PbiIndex(void);

    /// \}

public:
    /// \name PBI General Attributes
    /// \{

    /// \returns true if index has BarcodeData section
    bool HasBarcodeData(void) const;

    /// \returns true if index has MappedData section
    bool HasMappedData(void) const;

    /// \returns true if index has ReferenceData section
    bool HasReferenceData(void) const;

    /// \returns true if index has \b section
    /// \param[in] section PbiFile::Section identifier
    ///
    bool HasSection(const PbiFile::Section section) const;

    /// \returns index filename ("*.pbi")
    ///
    /// \note Returns an empty string if the underlying data was generated, not loaded from file.
    ///
    std::string Filename(void) const;

    /// \returns enum flags representing the file sections present
    PbiFile::Sections FileSections(void) const;

    /// \returns the number of records in the PBI (& associated BAM)
    uint32_t NumReads(void) const;

    /// \returns the PBI file's version number
    PbiFile::VersionEnum Version(void) const;

    /// \}

public:
    /// \name Lookup Data Components
    /// \{

    /// \returns barcode data lookup structure (may be empty, check result of HasBarcodeData)
    const BarcodeLookupData& BarcodeData(void) const;

    /// \returns basic data lookup structure
    const BasicLookupData& BasicData(void) const;

    /// \returns mapped data lookup structure (may be empty, check result of HasMappedData)
    const MappedLookupData& MappedData(void) const;

    /// \returns reference data lookup structure (may be empty, check result of HasReferenceData)
    const ReferenceLookupData& ReferenceData(void) const;

    /// }

private:
    PbiIndex(void);
    std::unique_ptr<internal::PbiIndexPrivate> d_;
};

} // namespace BAM
} // namespace PacBio

#include "internal/PbiIndex.inl"

#endif // PBIINDEX_H

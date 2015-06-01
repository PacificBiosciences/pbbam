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

#ifndef EXTERNALDATAREFERENCE_H
#define EXTERNALDATAREFERENCE_H

#include "pbbam/BamFile.h"
#include "pbbam/Config.h"
#include "pbbam/internal/DataSetElement.h"
#include "pbbam/internal/DataSetListElement.h"
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class PBBAM_EXPORT ExternalDataReference : public internal::DataSetElement
{
public:
    /// \name Constructors & Related Methods
    /// \{

    ExternalDataReference(void);

    /// Creates an external data reference pointing to \p BamFile
    ///
    /// Currently this means simply setting the MetaType to "SubreadFile.SubreadBamFile"
    /// and ResourceId to BamFile::Filename()
    ///
    ExternalDataReference(const BamFile& bamFile);

    using DataSetElement::DataSetElement;

    /// \}

public:
    /// \name Attributes
    ///

    /// \returns "Description" attribute value (or empty string if none present)
    ///
    const std::string& Description(void) const;

    /// \returns "Name" attribute value (or empty string if none present)
    ///
    const std::string& Name(void) const;

    /// \returns "MetaType" attribute value (or empty string if none present)
    ///
    const std::string& MetaType(void) const;

    /// \returns "ResourceId" attribute value (or empty string if none present)
    ///
    const std::string& ResourceId(void) const;

    /// \returns "Tags" attribute value (or empty string if none present)
    ///
    const std::string& Tags(void) const;

    /// \}

public:
    /// \name BAM Interoperability
    /// \{

    /// \returns true if reference points to a ".bam" file
    ///
    bool IsBamFile(void) const;

    /// \returns BamFile object representing this external data reference
    /// \exception throws if not a BAM file or fails to create BamFile object
    ///
    BamFile ToBamFile(void) const;

    /// \}

public:
    /// \name Attributes
    /// \{

    /// Sets "Description" attribute.
    ///
    /// \param[in] description
    /// \returns reference to this external data reference
    ///
    ExternalDataReference& Description(const std::string& description);

    /// Sets "MetaType" attribute.
    ///
    /// \param[in] metatype
    /// \returns reference to this external data reference
    ///
    ExternalDataReference& MetaType(const std::string& metatype);

    /// Sets "Name" attribute.
    ///
    /// \param[in] name
    /// \returns reference to this external data reference
    ///
    ExternalDataReference& Name(const std::string& name);

    /// Sets "ResourceId" attribute.
    ///
    /// \param[in] id
    /// \returns reference to this external data reference
    ///
    ExternalDataReference& ResourceId(const std::string& id);

    /// Sets "Tags" attribute.
    ///
    /// \param[in] tags
    /// \returns reference to this external data reference
    ///
    ExternalDataReference& Tags(const std::string& tags);

    /// \}
};

class PBBAM_EXPORT ExternalDataReferences : public internal::DataSetListElement<ExternalDataReference>
{
public:

    /// \name Constructors & Related Methods
    /// \{

    /// Creates empty external data reference list
    ExternalDataReferences(void);

    using DataSetListElement::DataSetListElement;

    /// \}

public:
    /// \n BAM-Specific Methods
    /// \{

    /// \returns list of BamFile objects for BAM files in this list
    std::vector<BamFile> BamFiles(void) const;

    /// \}

public:
    /// \name List Modification

    /// Adds \p ref to list
    ///
    /// \param[in] ref
    /// \returns reference to this list
    ///
    ExternalDataReferences& AddExternalRef(const ExternalDataReference& ref);

    /// Removes \p ref from list
    ///
    /// \param[in] ref
    /// \returns reference to this list
    ///
    ExternalDataReferences& RemoveExternalRef(const ExternalDataReference& ref);

    /// \}
};

} // namespace BAM
} // namespace PacBio

#endif // EXTERNALDATAREFERENCE_H

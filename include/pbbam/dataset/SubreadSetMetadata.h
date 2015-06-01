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

#ifndef SUBREADSETMETADATA_H
#define SUBREADSETMETADATA_H

#include "pbbam/dataset/DataSetMetadataBase.h"
#include "pbbam/internal/DataSetListElement.h"

namespace PacBio {
namespace BAM {

class BioSampleReferencesMetadata;
class BioSampleMetadata;
class BioSamplesMetadata;
class CollectionMetadata;
class CollectionsMetadata;
class CopyFilesMetadata;
class PrimaryMetadata;
class RunDetailsMetadata;
class SubreadSetMetadata;
class WellSampleMetadata;

class PBBAM_EXPORT BioSampleReferencesMetadata : public internal::DataSetElement
{
public:
    BioSampleReferencesMetadata(void);
    using DataSetElement::DataSetElement;
};

class PBBAM_EXPORT BioSampleMetadata : public internal::DataSetElement
{
public:
    BioSampleMetadata(void);
    using DataSetElement::DataSetElement;

public:
    const std::string& CreatedAt(void) const;
    const std::string& UniqueId(void) const;

public:
    BioSampleMetadata& CreatedAt(const std::string& createdAt);
    BioSampleMetadata& UniqueId(const std::string& uuid);
};

class PBBAM_EXPORT BioSamplesMetadata : public internal::DataSetListElement<BioSampleMetadata>
{
public:
    BioSamplesMetadata(void);
    using DataSetListElement::DataSetListElement;

public:
    BioSamplesMetadata& AddBioSample(const BioSampleMetadata& bioSample);
    BioSamplesMetadata& RemoveBioSample(const BioSampleMetadata& bioSample);
};

class PBBAM_EXPORT CollectionMetadata : public internal::DataSetElement
{
public:
    CollectionMetadata(void);
    using DataSetElement::DataSetElement;

public:
    const std::string& AutomationName(void) const;
    const std::string& CellIndex(void) const;
    const std::string& CellPac(void) const;
    const std::string& Context(void) const;
    const std::string& InstrCtrlVer(void) const;
    const std::string& InstrumentId(void) const;
    const std::string& InstrumentName(void) const;
    const std::string& SigProcVer(void) const;

    const PrimaryMetadata& Primary(void) const;
    const RunDetailsMetadata& RunDetails(void) const;
    const WellSampleMetadata& WellSample(void) const;

public:
    CollectionMetadata& AutomationName(const std::string& name);
    CollectionMetadata& CellIndex(const std::string& index);
    CollectionMetadata& CellPac(const std::string& pac);
    CollectionMetadata& Context(const std::string& context);
    CollectionMetadata& InstrCtrlVer(const std::string& ver);
    CollectionMetadata& InstrumentId(const std::string& id);
    CollectionMetadata& InstrumentName(const std::string& name);
    CollectionMetadata& SigProcVer(const std::string& ver);

    PrimaryMetadata& Primary(void);
    RunDetailsMetadata& RunDetails(void);
    WellSampleMetadata& WellSample(void);
};

class PBBAM_EXPORT CollectionsMetadata : public internal::DataSetListElement<CollectionMetadata>
{
public:
    CollectionsMetadata(void);
    using DataSetListElement::DataSetListElement;

public:
    CollectionsMetadata& AddCollection(const CollectionMetadata& collection);
    CollectionsMetadata& RemoveCollection(const CollectionMetadata& collection);
};

class PBBAM_EXPORT CopyFilesMetadata : public internal::DataSetElement
{
public:
    CopyFilesMetadata(void);
    using DataSetElement::DataSetElement;
};

class PBBAM_EXPORT PrimaryMetadata : public internal::DataSetElement
{
public:
    PrimaryMetadata(void);
    using DataSetElement::DataSetElement;

public:
    const std::string& AutomationName(void) const;
    const std::string& CollectionPathUri(void) const;
    const std::string& ContigFileName(void) const;
    const std::string& ResultsFolder(void) const;
    const std::string& SequencingCondition(void) const;

    const CopyFilesMetadata& CopyFiles(void) const;

public:
    PrimaryMetadata& AutomationName(const std::string& name);
    PrimaryMetadata& CollectionPathUri(const std::string& uri);
    PrimaryMetadata& ContigFileName(const std::string& name);
    PrimaryMetadata& ResultsFolder(const std::string& folder);
    PrimaryMetadata& SequencingCondition(const std::string& condition);

    CopyFilesMetadata& CopyFiles(void);
};

class PBBAM_EXPORT RunDetailsMetadata : public internal::DataSetElement
{
public:
    RunDetailsMetadata(void);
    using DataSetElement::DataSetElement;
public:
    const std::string& Name(void) const;
    const std::string& RunId(void) const;
public:
    RunDetailsMetadata& Name(const std::string& name);
    RunDetailsMetadata& RunId(const std::string& runId);
};

class PBBAM_EXPORT SubreadSetMetadata : public DataSetMetadataBase
{
public:
    SubreadSetMetadata(void);
    using DataSetMetadataBase::DataSetMetadataBase;

public:
    const BioSamplesMetadata& BioSamples(void) const;
    const CollectionsMetadata& Collections(void) const;

public:
    BioSamplesMetadata& BioSamples(void);
    CollectionsMetadata& Collections(void);
};

class PBBAM_EXPORT WellSampleMetadata : public internal::DataSetElement
{
public:
    WellSampleMetadata(void);
    using DataSetElement::DataSetElement;

public:
    const std::string& Comments(void) const;
    const std::string& Concentration(void) const;
    const std::string& PlateId(void) const;
    const std::string& SampleReuseEnabled(void) const;
    const std::string& SizeSelectionEnabled(void) const;
    const std::string& StageHotstartEnabled(void) const;
    const std::string& UniqueId(void) const;
    const std::string& UseCount(void) const;
    const std::string& WellName(void) const;

    const BioSampleReferencesMetadata& BioSampleReferences(void) const;

public:
    WellSampleMetadata& Comments(const std::string& comments);
    WellSampleMetadata& Concentration(const std::string& concentration);
    WellSampleMetadata& PlateId(const std::string& plateId);
    WellSampleMetadata& SampleReuseEnabled(const std::string& enabled);
    WellSampleMetadata& SizeSelectionEnabled(const std::string& enabled);
    WellSampleMetadata& StageHotstartEnabled(const std::string& enabled);
    WellSampleMetadata& UniqueId(const std::string& uuid);
    WellSampleMetadata& UseCount(const std::string& count);
    WellSampleMetadata& WellName(const std::string& name);

    BioSampleReferencesMetadata& BioSampleReferences(void);
};

} // namespace BAM
} // namespace PacBio

#endif // SUBREADSETMETADATA_H

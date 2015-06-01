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

#include "pbbam/dataset/ExternalDataReferences.h"
#include <boost/algorithm/string.hpp>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace PacBio::BAM::internal;
using namespace std;

// --------------------------------------
// ExternalDataReference implementation
// --------------------------------------

ExternalDataReference::ExternalDataReference(void)
    : DataSetElement("ExternalDataReference")
{ }

ExternalDataReference::ExternalDataReference(const BamFile& bamFile)
    : DataSetElement("ExternalDataReference")
{
    MetaType("SubreadFile.SubreadBamFile");
    ResourceId(bamFile.Filename());          // TODO: filename -> URI
}

const string &ExternalDataReference::Description(void) const
{ return Attribute("Description"); }

ExternalDataReference& ExternalDataReference::Description(const string& description)
{ Attribute("Description", description); return *this; }

bool ExternalDataReference::IsBamFile(void) const
{ return boost::algorithm::iends_with(ResourceId(), ".bam"); }

const string &ExternalDataReference::MetaType(void) const
{ return Attribute("MetaType"); }

ExternalDataReference& ExternalDataReference::MetaType(const string& metatype)
{ Attribute("MetaType", metatype); return *this; }

const string& ExternalDataReference::Name(void) const
{ return Attribute("Name"); }

ExternalDataReference& ExternalDataReference::Name(const string& name)
{ Attribute("Name", name); return *this; }

const string &ExternalDataReference::ResourceId(void) const
{ return Attribute("ResourceId"); }

ExternalDataReference& ExternalDataReference::ResourceId(const string& id)
{ Attribute("ResourceId", id); return *this; }

const string &ExternalDataReference::Tags(void) const
{ return Attribute("Tags"); }

ExternalDataReference& ExternalDataReference::Tags(const string& tags)
{ Attribute("Tags", tags); return *this; }

BamFile ExternalDataReference::ToBamFile(void) const
{
    if (!IsBamFile())
        throw std::exception();
    return BamFile(ResourceId()); // TODO: URI -> filename
}

// ---------------------------------------
// ExternalDataReferences implementation
// ---------------------------------------

ExternalDataReferences::ExternalDataReferences(void)
    : DataSetListElement<ExternalDataReference>("ExternalDataReferences")
{ }

ExternalDataReferences& ExternalDataReferences::AddExternalRef(const ExternalDataReference& ref)
{ AddChild(ref); return *this; }

ExternalDataReferences& ExternalDataReferences::RemoveExternalRef(const ExternalDataReference& ref)
{ RemoveChild(ref); return *this; }

vector<BamFile> ExternalDataReferences::BamFiles(void) const
{
    vector<BamFile> bamFiles;
    for ( const ExternalDataReference& ref : *this ) {
        if (ref.IsBamFile())
            bamFiles.push_back(ref.ToBamFile());
    }
    return bamFiles;
}

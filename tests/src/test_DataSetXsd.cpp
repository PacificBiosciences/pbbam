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

#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/DataSet.h>
#include <pbbam/DataSetXsd.h>

#include "PbbamTestData.h"

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(DataSetXsdTest, DefaultsOk)
{
    NamespaceRegistry registry;

    const NamespaceInfo& baseInfo = registry.Namespace(XsdType::BASE_DATA_MODEL);
    const NamespaceInfo& dsInfo   = registry.Namespace(XsdType::DATASETS);
    const NamespaceInfo& defaultInfo = registry.DefaultNamespace();

    EXPECT_EQ(XsdType::DATASETS, registry.DefaultXsd());

    EXPECT_EQ(string("pbds"),   dsInfo.Name());
    EXPECT_EQ(string("pbbase"), baseInfo.Name());
    EXPECT_EQ(string("pbds"),   defaultInfo.Name());

    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioBaseDataModel.xsd"), baseInfo.Uri());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"),      dsInfo.Uri());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"),      defaultInfo.Uri());
}

TEST(DataSetXsdTest, EditDefaultOk)
{
    NamespaceRegistry registry;
    registry.SetDefaultXsd(XsdType::DATASETS);

    const NamespaceInfo& defaultInfo = registry.DefaultNamespace();

    EXPECT_EQ(XsdType::DATASETS, registry.DefaultXsd());
    EXPECT_EQ(string("pbds"), defaultInfo.Name());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), defaultInfo.Uri());
}

TEST(DataSetXsdTest, EditRegistryOk)
{
    NamespaceRegistry registry;
    registry.Register(XsdType::DATASETS, NamespaceInfo("custom", "http://custom/uri.xsd"));

    const NamespaceInfo& dsInfo = registry.Namespace(XsdType::DATASETS);

    EXPECT_EQ(string("custom"),                dsInfo.Name());
    EXPECT_EQ(string("http://custom/uri.xsd"), dsInfo.Uri());
}

TEST(DataSetXsdTest, EditDatasetRegistry)
{
    DataSet dataset(DataSet::ALIGNMENT);
    dataset.CreatedAt("2015-01-27T09:00:01");
    dataset.MetaType("PacBio.DataSet.AlignmentSet");
    dataset.Name("DataSet_AlignmentSet");
    dataset.Tags("barcode moreTags mapping mytags");
    dataset.TimeStampedName("my_time_stamped_name");
    dataset.UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c");
    dataset.Attribute("xmlns",              "http://pacificbiosciences.com/PacBioDatasets.xsd")
           .Attribute("xmlns:xsi",          "http://www.w3.org/2001/XMLSchema-instance")
           .Attribute("xsi:schemaLocation", "http://pacificbiosciences.com/PacBioDatasets.xsd");

    ExternalResource ext("Fake.MetaType", "filename");
    ext.TimeStampedName("custom_tsn")
       .UniqueId("my_uuid");
    dataset.ExternalResources().Add(ext);

    dataset.Namespaces().Register(XsdType::BASE_DATA_MODEL, NamespaceInfo("custom", "http://custom/uri.xsd"));

    const string expectedXml =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:AlignmentSet CreatedAt=\"2015-01-27T09:00:01\" MetaType=\"PacBio.DataSet.AlignmentSet\" "
                "Name=\"DataSet_AlignmentSet\" Tags=\"barcode moreTags mapping mytags\" "
                "TimeStampedName=\"my_time_stamped_name\" "
                "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
                "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
                "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
                "xmlns:custom=\"http://custom/uri.xsd\" "
                "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t<custom:ExternalResources>\n"
        "\t\t<custom:ExternalResource MetaType=\"Fake.MetaType\" ResourceId=\"filename\" TimeStampedName=\"custom_tsn\" UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t</custom:ExternalResources>\n"
        "</pbds:AlignmentSet>\n";

    stringstream s;
    dataset.SaveToStream(s);
    EXPECT_EQ(expectedXml, s.str());
}

TEST(DataSetXsdTest, ElementRegistryOk)
{
    { // default namespaces

        DataSet ds;

        // append child elements that do not have a C++ built-in, nor namespace prefix with addition
        DataSetMetadata& metadata = ds.Metadata();
        metadata.AddChild(internal::DataSetElement("SummaryStats"));
        metadata.AddChild(internal::DataSetElement("CopyFiles"));
        metadata.AddChild(internal::DataSetElement("BioSamples"));
        metadata.AddChild(internal::DataSetElement("AutomationParameters"));

        stringstream s;
        ds.SaveToStream(s);
        const string output = s.str();

        // check that default namespace is propagated properly
        EXPECT_TRUE(output.find("pbds:SummaryStats") != string::npos);
        EXPECT_TRUE(output.find("pbmeta:CopyFiles") != string::npos);
        EXPECT_TRUE(output.find("pbsample:BioSamples") != string::npos);
        EXPECT_TRUE(output.find("pbbase:AutomationParameters") != string::npos);
    }

    { // custom namespaces

        DataSet ds;

        // setup custom namespaces
        ds.Namespaces().Register(XsdType::BASE_DATA_MODEL,     NamespaceInfo("custom_base",   "http://custom/base.xsd"));
        ds.Namespaces().Register(XsdType::COLLECTION_METADATA, NamespaceInfo("custom_meta",   "http://custom/meta.xsd"));
        ds.Namespaces().Register(XsdType::DATASETS,            NamespaceInfo("custom_ds",     "http://custom/datasets.xsd"));
        ds.Namespaces().Register(XsdType::SAMPLE_INFO,         NamespaceInfo("custom_sample", "http://custom/base.xsd"));

        // append child elements that do not have a C++ built-in, nor namespace prefix with addition
        DataSetMetadata& metadata = ds.Metadata();
        metadata.AddChild(internal::DataSetElement("SummaryStats"));
        metadata.AddChild(internal::DataSetElement("CopyFiles"));
        metadata.AddChild(internal::DataSetElement("BioSamples"));
        metadata.AddChild(internal::DataSetElement("AutomationParameters"));

        stringstream s;
        ds.SaveToStream(s);
        const string output = s.str();

        // check that custom namespace is propagated properly
        EXPECT_TRUE(output.find("custom_ds:SummaryStats") != string::npos);
        EXPECT_TRUE(output.find("custom_meta:CopyFiles") != string::npos);
        EXPECT_TRUE(output.find("custom_sample:BioSamples") != string::npos);
        EXPECT_TRUE(output.find("custom_base:AutomationParameters") != string::npos);
    }
}

// clang-format on

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

#ifdef PBBAM_TESTING
#define private public
#endif

#include "TestData.h"
#include <gtest/gtest.h>
#include <pbbam/DataSet.h>
#include <pbbam/internal/DataSetElement.h>
#include <stdexcept>
#include <sstream>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string ex2BamFn      = tests::Data_Dir + "/ex2.bam";
const string bamGroupFofn  = tests::Data_Dir + "/test_group_query/group.fofn";

const string ali1XmlFn = tests::Data_Dir + "/dataset/ali1.xml";
const string ali2XmlFn = tests::Data_Dir + "/dataset/ali2.xml";
const string ali3XmlFn = tests::Data_Dir + "/dataset/ali3.xml";
const string ali4XmlFn = tests::Data_Dir + "/dataset/ali4.xml";
const string mappingStaggeredXmlFn = tests::Data_Dir + "/dataset/bam_mapping_staggered.xml";
const string barcodeXmlFn = tests::Data_Dir + "/dataset/barcode.dataset.xml";
const string ccsReadXmlFn = tests::Data_Dir + "/dataset/ccsread.dataset.xml";
const string datasetFofn  = tests::Data_Dir + "/dataset/fofn.fofn";
const string hdfSubreadXmlFn    = tests::Data_Dir + "/dataset/hdfsubread_dataset.xml";
const string lambdaContigsXmlFn = tests::Data_Dir + "/dataset/lambda_contigs.xml";
const string pbalchemyXmlFn   = tests::Data_Dir + "/dataset/pbalchemy10kbp.xml";
const string referenceXmlFn   = tests::Data_Dir + "/dataset/reference.dataset.xml";
const string subread1XmlFn    = tests::Data_Dir + "/dataset/subread_dataset1.xml";
const string subread2XmlFn    = tests::Data_Dir + "/dataset/subread_dataset2.xml";
const string subread3XmlFn    = tests::Data_Dir + "/dataset/subread_dataset3.xml";
const string transformedXmlFn = tests::Data_Dir + "/dataset/transformed_rs_subread_dataset.xml";

static void TestAli1Xml(void);
static void TestAli2Xml(void);
static void TestAli3Xml(void);
static void TestAli4Xml(void);
static void TestMappingStaggeredXml(void);
static void TestBarcodeXml(void);
static void TestCcsReadXml(void);
static void TestHdfSubreadXml(void);
static void TestLambdaContigsXml(void);
static void TestPbalchemyXml(void);
static void TestReferenceXml(void);
static void TestSubread1Xml(void);
static void TestSubread2Xml(void);
static void TestSubread3Xml(void);
static void TestTransformedXml(void);

TEST(DataSetIOTest, FromBamFilename)
{
    DataSet dataset(ex2BamFn);

    EXPECT_EQ(1, dataset.ExternalResources().Size());
    const ExternalResource& bamRef = dataset.ExternalResources()[0];

    EXPECT_EQ(ex2BamFn, bamRef.ResourceId());
}

TEST(DataSetIOTest, FromBamFileObject)
{
    BamFile bamFile(ex2BamFn);
    DataSet dataset(bamFile.Filename());

    EXPECT_EQ(1, dataset.ExternalResources().Size());
    const ExternalResource& bamRef = dataset.ExternalResources()[0];

    EXPECT_EQ(ex2BamFn, bamRef.ResourceId());
}

TEST(DataSetIOTest, FromFofn)
{
    DataSet dataset(bamGroupFofn);
    EXPECT_EQ(3, dataset.ExternalResources().Size());
}

TEST(DataSetIOTest, FromXml)
{
    EXPECT_NO_THROW(TestAli1Xml());
    EXPECT_NO_THROW(TestAli2Xml());
    EXPECT_NO_THROW(TestAli3Xml());
    EXPECT_NO_THROW(TestAli4Xml());
    EXPECT_NO_THROW(TestMappingStaggeredXml());
    EXPECT_NO_THROW(TestBarcodeXml());
    EXPECT_NO_THROW(TestCcsReadXml());
    EXPECT_NO_THROW(TestHdfSubreadXml());
    EXPECT_NO_THROW(TestLambdaContigsXml());
    EXPECT_NO_THROW(TestPbalchemyXml());
    EXPECT_NO_THROW(TestReferenceXml());
    EXPECT_NO_THROW(TestSubread1Xml());
    EXPECT_NO_THROW(TestSubread2Xml());
    EXPECT_NO_THROW(TestSubread3Xml());
    EXPECT_NO_THROW(TestTransformedXml());
}

TEST(DataSetIOTest, ThrowsOnNonexistentFofnFile)
{
    EXPECT_THROW(
    {
        DataSet dataset("does/not/exist.fofn");

    }, std::exception);
}

TEST(DataSetIOTest, ThrowsOnNonexistentXmlFile)
{
    EXPECT_THROW(
    {
        DataSet dataset("does/not/exist.xml");

    }, std::exception);
}

TEST(DataSetIOTest, ToXml)
{
    // top-level data
    DataSet dataset(DataSet::ALIGNMENT);
    dataset.CreatedAt("2015-01-27T09:00:01")
           .MetaType("PacBio.DataSet.AlignmentSet")
           .Name("DataSet_AlignmentSet")
           .Tags("barcode moreTags mapping mytags")
           .UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c")
           .Version("2.3.0");
    dataset.Attribute("xmlns","http://pacificbiosciences.com/PacBioDataModel.xsd")
           .Attribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
           .Attribute("xsi:schemaLocation",
        "http://pacificbiosciences.com/PacBioDataModel.xsd "
        "../../../../../../../../common/datamodel/SequEl/EndToEnd/xsd/PacBioSecondaryDataModel.xsd");

    // external resources
    ExternalResource resource1;
    resource1.Name("Third Alignments BAM")
        .Description("Points to an example Alignments BAM file.")
        .MetaType("AlignmentFile.AlignmentBamFile")
        .ResourceId("file:/mnt/path/to/alignments2.bam")
        .Tags("Example");
    FileIndex pbi1;
    pbi1.MetaType("PacBio.Index.PacBioIndex");
    pbi1.ResourceId("file:/mnt/path/to/alignments2.pbi");
    resource1.FileIndices().Add(pbi1);
    dataset.ExternalResources().Add(resource1);

    ExternalResource resource2;
    resource2.Name("Fourth Alignments BAM")
        .Description("Points to another example Alignments BAM file, by relative path.")
        .MetaType("AlignmentFile.AlignmentBamFile")
        .ResourceId("file:./alignments3.bam")
        .Tags("Example");
    FileIndex pbi2;
    pbi2.MetaType("PacBio.Index.PacBioIndex");
    pbi2.ResourceId("file:/mnt/path/to/alignments3.pbi");
    resource2.FileIndices().Add(pbi2);
    dataset.ExternalResources().Add(resource2);

    // sub-datasets with filters
    DataSetBase subDataSet1;
    subDataSet1.Name("HighQuality Read Alignments")
               .UniqueId("ab95d0a3-94b8-4918-b3af-a3f81bbe519c")
               .Version("2.3.0");
    Filter filter1;
    filter1.Properties().Add(Property("rq", "0.85", ">"));
    subDataSet1.Filters().Add(filter1);
    dataset.SubDataSets().Add(subDataSet1);

    DataSetBase subDataSet2;
    subDataSet2.Name("Alignments to chromosome 1")
               .UniqueId("ac95d0a3-94b8-4918-b3af-a3f81bbe519c")
               .Version("2.3.0");
    Filter filter2;
    filter2.Properties().Add(Property("RNAME", "chr1", "=="));
    subDataSet2.Filters().Add(filter2);
    dataset.SubDataSets().Add(subDataSet2);

    // write dataset
    const string expectedXml =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<AlignmentSet CreatedAt=\"2015-01-27T09:00:01\" MetaType=\"PacBio.DataSet.AlignmentSet\" "
                "Name=\"DataSet_AlignmentSet\" Tags=\"barcode moreTags mapping mytags\" "
                "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"2.3.0\" "
                "xmlns=\"http://pacificbiosciences.com/PacBioDataModel.xsd\" "
                "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDataModel.xsd "
                "../../../../../../../../common/datamodel/SequEl/EndToEnd/xsd/PacBioSecondaryDataModel.xsd\">\n"
        "\t<ExternalResources>\n"
        "\t\t<ExternalResource Description=\"Points to an example Alignments BAM file.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" Name=\"Third Alignments BAM\" "
                "ResourceId=\"file:/mnt/path/to/alignments2.bam\" Tags=\"Example\">\n"
        "\t\t\t<FileIndices>\n"
        "\t\t\t\t<FileIndex MetaType=\"PacBio.Index.PacBioIndex\" ResourceId=\"file:/mnt/path/to/alignments2.pbi\" />\n"
        "\t\t\t</FileIndices>\n"
        "\t\t</ExternalResource>\n"
        "\t\t<ExternalResource Description=\"Points to another example Alignments BAM file, by relative path.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" Name=\"Fourth Alignments BAM\" "
                "ResourceId=\"file:./alignments3.bam\" Tags=\"Example\">\n"
        "\t\t\t<FileIndices>\n"
        "\t\t\t\t<FileIndex MetaType=\"PacBio.Index.PacBioIndex\" ResourceId=\"file:/mnt/path/to/alignments3.pbi\" />\n"
        "\t\t\t</FileIndices>\n"
        "\t\t</ExternalResource>\n"
        "\t</ExternalResources>\n"
        "\t<DataSets>\n"
        "\t\t<DataSet Name=\"HighQuality Read Alignments\" UniqueId=\"ab95d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"2.3.0\">\n"
        "\t\t\t<Filters>\n"
        "\t\t\t\t<Filter>\n"
        "\t\t\t\t\t<Properties>\n"
        "\t\t\t\t\t\t<Property Name=\"rq\" Operator=\">\" Value=\"0.85\" />\n"
        "\t\t\t\t\t</Properties>\n"
        "\t\t\t\t</Filter>\n"
        "\t\t\t</Filters>\n"
        "\t\t</DataSet>\n"
        "\t\t<DataSet Name=\"Alignments to chromosome 1\" UniqueId=\"ac95d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"2.3.0\">\n"
        "\t\t\t<Filters>\n"
        "\t\t\t\t<Filter>\n"
        "\t\t\t\t\t<Properties>\n"
        "\t\t\t\t\t\t<Property Name=\"RNAME\" Operator=\"==\" Value=\"chr1\" />\n"
        "\t\t\t\t\t</Properties>\n"
        "\t\t\t\t</Filter>\n"
        "\t\t\t</Filters>\n"
        "\t\t</DataSet>\n"
        "\t</DataSets>\n"
        "</AlignmentSet>\n";

    stringstream s;
    dataset.SaveToStream(s);
    EXPECT_EQ(expectedXml, s.str());
}

static void TestAli1Xml(void)
{
    const DataSet dataset(ali1XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),                  dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.AlignmentSet"),          dataset.MetaType());
    EXPECT_EQ(string("DataSet_AlignmentSet"),                 dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"),      dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Alignments BAM file."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/alignments0.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments0.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Alignments BAM file, by relative path."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:./alignments1.bam"),  resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments1.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.85"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        }
        else {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("RNAME"), property.Name());
            EXPECT_EQ(string("chr1"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }
}

static void TestAli2Xml(void)
{
    const DataSet dataset(ali2XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),                  dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.AlignmentSet"),          dataset.MetaType());
    EXPECT_EQ(string("DataSet_AlignmentSet"),                 dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"),      dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Alignments BAM file."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/alignments2.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments2.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Alignments BAM file, by relative path."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:./alignments3.bam"),  resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments3.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.85"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        }
        else {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("RNAME"), property.Name());
            EXPECT_EQ(string("chr1"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }
}

static void TestAli3Xml(void)
{
    const DataSet dataset(ali3XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),                  dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.AlignmentSet"),          dataset.MetaType());
    EXPECT_EQ(string("DataSet_AlignmentSet"),                 dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"),      dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Alignments BAM file."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/alignments2.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments2.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Alignments BAM file, by relative path."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:./alignments3.bam"),  resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments3.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.75"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        }
        else {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("RNAME"), property.Name());
            EXPECT_EQ(string("chr1"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }
}

static void TestAli4Xml(void)
{
    const DataSet dataset(ali4XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),                  dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.AlignmentSet"),          dataset.MetaType());
    EXPECT_EQ(string("DataSet_AlignmentSet"),                 dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"),      dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Alignments BAM file."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/alignments0.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments0.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Alignments BAM file, by relative path."), resource.Description());
            EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:./alignments1.bam"),  resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/alignments1.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.85"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        }
        else {
            EXPECT_EQ(string(""), subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("RNAME"), property.Name());
            EXPECT_EQ(string("chr1"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }
}

static void TestMappingStaggeredXml(void)
{
    const DataSet dataset(mappingStaggeredXmlFn);
    EXPECT_EQ(DataSet::GENERIC, dataset.Type());
    EXPECT_EQ(string("2015-05-13T10:58:26"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.DataSet"), dataset.MetaType());
    EXPECT_EQ(string(""), dataset.Name());
    EXPECT_EQ(string(""), dataset.Tags());
    EXPECT_EQ(string("30f72098-bc5b-e06b-566c-8b28dda909a8"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string(""), resource.Name());
            EXPECT_EQ(string(""), resource.Description());
            EXPECT_EQ(string(""), resource.MetaType());
            EXPECT_EQ(string("file:tests/data/bam_mapping_1.bam"), resource.ResourceId());
            EXPECT_EQ(string(""), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:tests/data/bam_mapping_1.bam.bai"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string(""), resource.Name());
            EXPECT_EQ(string(""), resource.Description());
            EXPECT_EQ(string(""), resource.MetaType());
            EXPECT_EQ(string("file:tests/data/bam_mapping_2.bam"), resource.ResourceId());
            EXPECT_EQ(string(""), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:tests/data/bam_mapping_2.bam.bai"), index.ResourceId());
        }
    }

    // ?
    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(string("2015-05-13T10:58:26"),    subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string(""), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("c5402d06-4643-057c-e300-fe229b4e8909"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const ExternalResources& resources = subdataset.ExternalResources();
            ASSERT_EQ(1, resources.Size());
            const ExternalResource& resource = resources[0];
            EXPECT_EQ(string("file:tests/data/bam_mapping_2.bam"), resource.ResourceId());
            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:tests/data/bam_mapping_2.bam.bai"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("2015-05-13T10:58:26"),    subdataset.CreatedAt());
            EXPECT_EQ(string(""), subdataset.MetaType());
            EXPECT_EQ(string(""), subdataset.Name());
            EXPECT_EQ(string(""), subdataset.Tags());
            EXPECT_EQ(string("f8b54a55-5fb7-706f-ab35-39afc9c86924"), subdataset.UniqueId());
            EXPECT_EQ(string("2.3.0"), subdataset.Version());

            const ExternalResources& resources = subdataset.ExternalResources();
            ASSERT_EQ(1, resources.Size());
            const ExternalResource& resource = resources[0];
            EXPECT_EQ(string("file:tests/data/bam_mapping_1.bam"), resource.ResourceId());
            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:tests/data/bam_mapping_1.bam.bai"), index.ResourceId());
        }
    }
}

static void TestBarcodeXml(void)
{
    const DataSet dataset(barcodeXmlFn);
    EXPECT_EQ(DataSet::BARCODE, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.BarcodeSet"), dataset.MetaType());
    EXPECT_EQ(string("DataSet_BarcodeSet"), dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(string("First Barcodes FASTA"), resource.Name());
    EXPECT_EQ(string("Points to an example Barcodes FASTA file."), resource.Description());
    EXPECT_EQ(string("BarcodeFile.BarcodeFastaFile"), resource.MetaType());
    EXPECT_EQ(string("file:///mnt/path/to/barcode.fasta"), resource.ResourceId());
    EXPECT_EQ(string("Example"), resource.Tags());

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("30"),     metadata.NumRecords());
    EXPECT_EQ(string("400"),    metadata.TotalLength());

    // access metadata extensions directly for now
    EXPECT_EQ(string("paired"), metadata.ChildText("BarcodeConstruction"));
}

static void TestCcsReadXml(void)
{
    const DataSet dataset(ccsReadXmlFn);
    EXPECT_EQ(DataSet::CONSENSUS_READ, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.ConsensusReadSet"), dataset.MetaType());
    EXPECT_EQ(string("DataSet_ConsensusReadSet"), dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First ConsensusRead BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example ConsensusRead BAM file."), resource.Description());
            EXPECT_EQ(string("PacBio.ConsensusReadFile.ConsensusReadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/ccsreads0.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("PacBio.Index.PacBioIndex"), index.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/ccsreads0.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second ConsensusRead BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example ConsensusRead BAM file."), resource.Description());
            EXPECT_EQ(string("PacBio.ConsensusReadFile.ConsensusReadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/ccsreads1.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("PacBio.Index.PacBioIndex"), index.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/ccsreads0.pbi"), index.ResourceId());
        }
    }
}

static void TestHdfSubreadXml(void)
{
    // Looks like a bunch of TYPOS in XML file !!
}

static void TestLambdaContigsXml(void)
{
    const DataSet dataset(lambdaContigsXmlFn);
    EXPECT_EQ(DataSet::REFERENCE, dataset.Type());
    EXPECT_EQ(string("2015-05-28T10:56:36"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.ReferenceSet"), dataset.MetaType());
    EXPECT_EQ(string(""), dataset.Name());
    EXPECT_EQ(string(""), dataset.Tags());
    EXPECT_EQ(string("596e87db-34f9-d2fd-c905-b017543170e1"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(string("file:tests/data/lambda_contigs.fasta"), resource.ResourceId());
}

static void TestPbalchemyXml(void)
{
    const DataSet dataset(pbalchemyXmlFn);
    EXPECT_EQ(DataSet::GENERIC, dataset.Type());
    EXPECT_EQ(string("2015-05-22T16:56:16"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.DataSet"), dataset.MetaType());
    EXPECT_EQ(string(""), dataset.Name());
    EXPECT_EQ(string(""), dataset.Tags());
    EXPECT_EQ(string("58e3f7c5-24c1-b58b-fbd5-37de268cc2f0"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(string("file:tests/data/pbalchemy10kbp.pbalign.sorted.pbver1.bam"), resource.ResourceId());
    const FileIndices& fileIndices = resource.FileIndices();
    ASSERT_EQ(1, fileIndices.Size());
    const FileIndex& index = fileIndices[0];
    EXPECT_EQ(string("file:tests/data/pbalchemy10kbp.pbalign.sorted.pbver1.bam.bai"), index.ResourceId());

    // TYPOs: Should be Filter Properties/Property not Parameter(s)

}

static void TestReferenceXml(void)
{
    const DataSet dataset(referenceXmlFn);
    EXPECT_EQ(DataSet::REFERENCE, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.ReferenceSet"), dataset.MetaType());
    EXPECT_EQ(string("DataSet_ReferenceSet"), dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(string("First References FASTA"), resource.Name());
    EXPECT_EQ(string("Points to an example references FASTA file."), resource.Description());
    EXPECT_EQ(string("PacBio.ReferenceFile.ReferenceFastaFile"), resource.MetaType());
    EXPECT_EQ(string("file:///mnt/path/to/reference.fasta"), resource.ResourceId());
    EXPECT_EQ(string("Example"), resource.Tags());
    const FileIndices& fileIndices = resource.FileIndices();
    ASSERT_EQ(2, fileIndices.Size());
    for (size_t i = 0; i < fileIndices.Size(); ++i) {
        const FileIndex& index = fileIndices[i];
        if (i == 0) {
            EXPECT_EQ(string("PacBio.Index.SaWriterIndex"), index.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/reference.fasta.sa"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("PacBio.Index.SamIndex"), index.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/reference.fasta.fai"), index.ResourceId());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("500"),     metadata.NumRecords());
    EXPECT_EQ(string("5000000"), metadata.TotalLength());

    // access metadata extensions directly for now
    EXPECT_EQ(string("Tribble"), metadata.ChildText("Organism"));
    EXPECT_EQ(string("Diploid"), metadata.ChildText("Ploidy"));

    const internal::DataSetListElement<internal::DataSetElement>& contigs =
            metadata.Child<internal::DataSetListElement<internal::DataSetElement> >("Contigs");
    ASSERT_EQ(1, contigs.NumChildren());
    const internal::DataSetElement& contig = contigs[0];
    EXPECT_EQ(string("gi|229359445|emb|AM181176.4|"), contig.Attribute("Name"));
    EXPECT_EQ(string("Pseudomonas fluorescens SBW25 complete genome|quiver"), contig.Attribute("Description"));
    EXPECT_EQ(string("6722109"), contig.Attribute("Length"));
    EXPECT_EQ(string("f627c795efad7ce0050ed42b942d408e"), contig.Attribute("Digest"));
}

static void TestSubread1Xml(void)
{
    const DataSet dataset(subread1XmlFn);
    EXPECT_EQ(DataSet::SUBREAD, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(string("DataSet_SubreadSet"), dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Subreads BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Subreads BAM file."), resource.Description());
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/subreads0.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/subreads0.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Subreads BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Subreads BAM file."), resource.Description());
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/subreads1.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/subreads0.pbi"), index.ResourceId());
        }
    }

    const Filters& filters = dataset.Filters();
    ASSERT_EQ(2, filters.Size());
    for (size_t i = 0; i < filters.Size(); ++i) {
        const Filter& filter = filters[i];
        if (i == 0) {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.75"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        } else {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("QNAME"), property.Name());
            EXPECT_EQ(string("100/0/0_100"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("500"),    metadata.NumRecords());
    EXPECT_EQ(string("500000"), metadata.TotalLength());
}

static void TestSubread2Xml(void)
{
    const DataSet dataset(subread2XmlFn);
    EXPECT_EQ(DataSet::SUBREAD, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),    dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(string("DataSet_SubreadSet"), dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Subreads BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Subreads BAM file."), resource.Description());
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/subreads2.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/subreads2.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Subreads BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Subreads BAM file."), resource.Description());
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/subreads3.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/subreads3.pbi"), index.ResourceId());
        }
    }

    const Filters& filters = dataset.Filters();
    ASSERT_EQ(2, filters.Size());
    for (size_t i = 0; i < filters.Size(); ++i) {
        const Filter& filter = filters[i];
        if (i == 0) {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.75"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        } else {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("QNAME"), property.Name());
            EXPECT_EQ(string("100/0/0_100"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("500"),    metadata.NumRecords());
    EXPECT_EQ(string("500000"), metadata.TotalLength());
}

static void TestSubread3Xml(void)
{
    const DataSet dataset(subread3XmlFn);
    EXPECT_EQ(DataSet::SUBREAD, dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(string("DataSet_SubreadSet"), dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),         dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xsi:schemaLocation"));

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("First Subreads BAM"), resource.Name());
            EXPECT_EQ(string("Points to an example Subreads BAM file."), resource.Description());
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/subreads2.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/subreads2.pbi"), index.ResourceId());
        }
        else {
            EXPECT_EQ(string("Second Subreads BAM"), resource.Name());
            EXPECT_EQ(string("Points to another example Subreads BAM file."), resource.Description());
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/path/to/subreads3.bam"), resource.ResourceId());
            EXPECT_EQ(string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(string("file:///mnt/path/to/subreads3.pbi"), index.ResourceId());
        }
    }

    const Filters& filters = dataset.Filters();
    ASSERT_EQ(2, filters.Size());
    for (size_t i = 0; i < filters.Size(); ++i) {
        const Filter& filter = filters[i];
        if (i == 0) {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("rq"), property.Name());
            EXPECT_EQ(string("0.85"), property.Value());
            EXPECT_EQ(string(">"), property.Operator());
        } else {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(string("QNAME"), property.Name());
            EXPECT_EQ(string("100/0/0_100"), property.Value());
            EXPECT_EQ(string("=="), property.Operator());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("500"),    metadata.NumRecords());
    EXPECT_EQ(string("500000"), metadata.TotalLength());
}

static void TestTransformedXml(void)
{
    const DataSet dataset(transformedXmlFn);
    EXPECT_EQ(DataSet::HDF_SUBREAD, dataset.Type());
    EXPECT_EQ(string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(string("Subreads from run r001173_42129_130607"), dataset.Name());
    EXPECT_EQ(string("pacbio.secondary.instrument=RS"), dataset.Tags());
    EXPECT_EQ(string("abbc9183-b01e-4671-8c12-19efee534647"), dataset.UniqueId());
    EXPECT_EQ(string("0.5"), dataset.Version());
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDataModel.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema"),         dataset.Attribute("xmlns:xs"));
    EXPECT_EQ(string("http://www.w3.org/2005/xpath-functions"), dataset.Attribute("xmlns:fn"));
    EXPECT_EQ(string("java:java.util.UUID"), dataset.Attribute("xmlns:uuid"));
    EXPECT_EQ(string("http://whatever"), dataset.Attribute("xmlns:bax"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(3, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(string("PacBio.SubreadFile.BaxFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/secondary-siv/testdata/LIMS/2590727/0001/Analysis_Results/m130608_033634_42129_c100515232550000001823076608221351_s1_p0.0.bax.h5"),
                      resource.ResourceId());
        }
        else if (i == 1) {
            EXPECT_EQ(string("PacBio.SubreadFile.BaxFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/secondary-siv/testdata/LIMS/2590727/0001/Analysis_Results/m130608_033634_42129_c100515232550000001823076608221351_s1_p0.1.bax.h5"),
                      resource.ResourceId());
        }
        else {
            EXPECT_EQ(string("PacBio.SubreadFile.BaxFile"), resource.MetaType());
            EXPECT_EQ(string("file:///mnt/secondary-siv/testdata/LIMS/2590727/0001/Analysis_Results/m130608_033634_42129_c100515232550000001823076608221351_s1_p0.2.bax.h5"),
                      resource.ResourceId());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("150000"),   metadata.NumRecords());
    EXPECT_EQ(string("50000000"), metadata.TotalLength());
}


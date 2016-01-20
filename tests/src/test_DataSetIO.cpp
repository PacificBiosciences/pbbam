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
#include "../src/FileUtils.h"
#include <gtest/gtest.h>
#include <pbbam/DataSet.h>
#include <pbbam/internal/DataSetElement.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
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

static void TestFromXmlString(void);
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

static inline
void changeCurrentDirectory(const std::string& dir)
{ ASSERT_EQ(0, chdir(dir.c_str())); }

TEST(DataSetIOTest, FromBamFilename)
{
    DataSet dataset(ex2BamFn);

    EXPECT_EQ(1, dataset.ExternalResources().Size());
    const ExternalResource& bamRef = dataset.ExternalResources()[0];

    EXPECT_EQ(ex2BamFn, bamRef.ResourceId());
}

TEST(DataSetIOTest, FromBamFilenames)
{
    std::ifstream fofn(bamGroupFofn);
    std::vector<std::string> files;
    std::string file;
    while (std::getline(fofn, file)) if (!file.empty()) files.emplace_back(file);
    DataSet dataset(files);
    EXPECT_EQ(3, dataset.ExternalResources().Size());
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
    EXPECT_NO_THROW(TestFromXmlString());
}

TEST(DataSetIOTest, FromXmlFile)
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
    dataset.CreatedAt("2015-01-27T09:00:01");
    dataset.MetaType("PacBio.DataSet.AlignmentSet");
    dataset.Name("DataSet_AlignmentSet");
    dataset.Tags("barcode moreTags mapping mytags");
    dataset.TimeStampedName("my_tsn");
    dataset.UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c");
    dataset.Attribute("xmlns",              "http://pacificbiosciences.com/PacBioDatasets.xsd")
           .Attribute("xmlns:xsi",          "http://www.w3.org/2001/XMLSchema-instance")
           .Attribute("xsi:schemaLocation", "http://pacificbiosciences.com/PacBioDatasets.xsd");

    // external resources
    ExternalResource resource1("AlignmentFile.AlignmentBamFile", "file:/mnt/path/to/alignments2.bam");    
    resource1.Name("Third Alignments BAM");
    resource1.Description("Points to an example Alignments BAM file.");
    resource1.Tags("Example");
    resource1.TimeStampedName("my_tsn");
    resource1.UniqueId("my_uuid");
    FileIndex pbi1("PacBio.Index.PacBioIndex", "file:/mnt/path/to/alignments2.pbi");
    pbi1.TimeStampedName("my_tsn");
    pbi1.UniqueId("my_uuid");
    resource1.FileIndices().Add(pbi1);
    dataset.ExternalResources().Add(resource1);

    ExternalResource resource2("AlignmentFile.AlignmentBamFile", "file:./alignments3.bam");
    resource2.Name("Fourth Alignments BAM");
    resource2.Description("Points to another example Alignments BAM file, by relative path.");
    resource2.Tags("Example");
    resource2.TimeStampedName("my_tsn");
    resource2.UniqueId("my_uuid");
    FileIndex pbi2("PacBio.Index.PacBioIndex", "file:/mnt/path/to/alignments3.pbi");
    pbi2.TimeStampedName("my_tsn");
    pbi2.UniqueId("my_uuid");

    resource2.FileIndices().Add(pbi2);
    dataset.ExternalResources().Add(resource2);

    // sub-datasets with filters
    DataSetBase subDataSet1;
    subDataSet1.Name("HighQuality Read Alignments");
    subDataSet1.TimeStampedName("my_tsn");
    subDataSet1.UniqueId("ab95d0a3-94b8-4918-b3af-a3f81bbe519c");
    Filter filter1;
    filter1.Properties().Add(Property("rq", "0.85", ">"));
    subDataSet1.Filters().Add(filter1);
    dataset.SubDataSets().Add(subDataSet1);

    DataSetBase subDataSet2;
    subDataSet2.Name("Alignments to chromosome 1");
    subDataSet2.TimeStampedName("my_tsn");
    subDataSet2.UniqueId("ac95d0a3-94b8-4918-b3af-a3f81bbe519c");
    Filter filter2;
    filter2.Properties().Add(Property("RNAME", "chr1", "=="));
    subDataSet2.Filters().Add(filter2);
    dataset.SubDataSets().Add(subDataSet2);

    // write dataset
    const string expectedXml =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:AlignmentSet "
                "CreatedAt=\"2015-01-27T09:00:01\" "
                "MetaType=\"PacBio.DataSet.AlignmentSet\" "
                "Name=\"DataSet_AlignmentSet\" "
                "Tags=\"barcode moreTags mapping mytags\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
                "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
                "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
                "xmlns:pbbase=\"http://pacificbiosciences.com/PacBioBaseDataModel.xsd\" "
                "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t<pbbase:ExternalResources>\n"
        "\t\t<pbbase:ExternalResource "
                "Description=\"Points to an example Alignments BAM file.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" "
                "Name=\"Third Alignments BAM\" "
                "ResourceId=\"file:/mnt/path/to/alignments2.bam\" "
                "Tags=\"Example\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"my_uuid\" Version=\"3.0.1\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
                "MetaType=\"PacBio.Index.PacBioIndex\" "
                "ResourceId=\"file:/mnt/path/to/alignments2.pbi\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t\t<pbbase:ExternalResource "
                "Description=\"Points to another example Alignments BAM file, by relative path.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" "
                "Name=\"Fourth Alignments BAM\" "
                "ResourceId=\"file:./alignments3.bam\" "
                "Tags=\"Example\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"my_uuid\" Version=\"3.0.1\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
                "MetaType=\"PacBio.Index.PacBioIndex\" "
                "ResourceId=\"file:/mnt/path/to/alignments3.pbi\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t</pbbase:ExternalResources>\n"
        "\t<pbds:DataSets>\n"
        "\t\t<pbds:DataSet "
                "MetaType=\"PacBio.DataSet.DataSet\" "
                "Name=\"HighQuality Read Alignments\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"ab95d0a3-94b8-4918-b3af-a3f81bbe519c\" "
                "Version=\"3.0.1\">\n"
        "\t\t\t<pbds:Filters>\n"
        "\t\t\t\t<pbds:Filter>\n"
        "\t\t\t\t\t<pbbase:Properties>\n"
        "\t\t\t\t\t\t<pbbase:Property Name=\"rq\" Operator=\">\" Value=\"0.85\" />\n"
        "\t\t\t\t\t</pbbase:Properties>\n"
        "\t\t\t\t</pbds:Filter>\n"
        "\t\t\t</pbds:Filters>\n"
        "\t\t</pbds:DataSet>\n"
        "\t\t<pbds:DataSet "
                "MetaType=\"PacBio.DataSet.DataSet\" "
                "Name=\"Alignments to chromosome 1\" "
                "TimeStampedName=\"my_tsn\" "
                "UniqueId=\"ac95d0a3-94b8-4918-b3af-a3f81bbe519c\" "
                "Version=\"3.0.1\">\n"
        "\t\t\t<pbds:Filters>\n"
        "\t\t\t\t<pbds:Filter>\n"
        "\t\t\t\t\t<pbbase:Properties>\n"
        "\t\t\t\t\t\t<pbbase:Property Name=\"RNAME\" Operator=\"==\" Value=\"chr1\" />\n"
        "\t\t\t\t\t</pbbase:Properties>\n"
        "\t\t\t\t</pbds:Filter>\n"
        "\t\t\t</pbds:Filters>\n"
        "\t\t</pbds:DataSet>\n"
        "\t</pbds:DataSets>\n"
        "</pbds:AlignmentSet>\n";

    stringstream s;
    dataset.SaveToStream(s);
    EXPECT_EQ(expectedXml, s.str());
}

static void TestFromXmlString(void)
{
    const string inputXml =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:AlignmentSet "
            "CreatedAt=\"2015-01-27T09:00:01\" "
            "MetaType=\"PacBio.DataSet.AlignmentSet\" "
            "Name=\"DataSet_AlignmentSet\" "
            "Tags=\"barcode moreTags mapping mytags\" "
            "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" "
            "Version=\"2.3.0\" "
            "xmlns=\"http://pacificbiosciences.com/PacBioDataModel.xsd\" "
            "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
            "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDataModel.xsd\">\n"
        "\t<pbbase:ExternalResources>\n"
        "\t\t<pbbase:ExternalResource "
                "Description=\"Points to an example Alignments BAM file.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" "
                "Name=\"Third Alignments BAM\" "
                "ResourceId=\"file:/mnt/path/to/alignments2.bam\" "
                "Tags=\"Example\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
                    "MetaType=\"PacBio.Index.PacBioIndex\" "
                    "ResourceId=\"file:/mnt/path/to/alignments2.pbi\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t\t<pbbase:ExternalResource "
                "Description=\"Points to another example Alignments BAM file, by relative path.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" "
                "Name=\"Fourth Alignments BAM\" "
                "ResourceId=\"file:./alignments3.bam\" "
                "Tags=\"Example\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
                    "MetaType=\"PacBio.Index.PacBioIndex\" "
                    "ResourceId=\"file:/mnt/path/to/alignments3.pbi\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t</pbbase:ExternalResources>\n"
        "\t<pbds:DataSets>\n"
        "\t\t<pbds:DataSet "
                "Name=\"HighQuality Read Alignments\" "
                "UniqueId=\"ab95d0a3-94b8-4918-b3af-a3f81bbe519c\" "
                "Version=\"2.3.0\">\n"
        "\t\t\t<pbds:Filters>\n"
        "\t\t\t\t<pbds:Filter>\n"
        "\t\t\t\t\t<pbbase:Properties>\n"
        "\t\t\t\t\t\t<pbbase:Property Name=\"rq\" Operator=\">\" Value=\"0.85\" />\n"
        "\t\t\t\t\t</pbbase:Properties>\n"
        "\t\t\t\t</pbds:Filter>\n"
        "\t\t\t</pbds:Filters>\n"
        "\t\t</pbds:DataSet>\n"
        "\t\t<pbds:DataSet "
                "Name=\"Alignments to chromosome 1\" "
                "UniqueId=\"ac95d0a3-94b8-4918-b3af-a3f81bbe519c\" "
                "Version=\"2.3.0\">\n"
        "\t\t\t<pbds:Filters>\n"
        "\t\t\t\t<pbds:Filter>\n"
        "\t\t\t\t\t<pbbase:Properties>\n"
        "\t\t\t\t\t\t<pbbase:Property Name=\"RNAME\" Operator=\"==\" Value=\"chr1\" />\n"
        "\t\t\t\t\t</pbbase:Properties>\n"
        "\t\t\t\t</pbds:Filter>\n"
        "\t\t\t</pbds:Filters>\n"
        "\t\t</pbds:DataSet>\n"
        "\t</pbds:DataSets>\n"
        "</pbds:AlignmentSet>\n";

    const DataSet dataset = DataSet::FromXml(inputXml);

    EXPECT_EQ(DataSet::ALIGNMENT,                     dataset.Type());
    EXPECT_EQ("2015-01-27T09:00:01",                  dataset.CreatedAt());
    EXPECT_EQ("PacBio.DataSet.AlignmentSet",          dataset.MetaType());
    EXPECT_EQ("DataSet_AlignmentSet",                 dataset.Name());
    EXPECT_EQ("barcode moreTags mapping mytags",      dataset.Tags());
    EXPECT_EQ("b095d0a3-94b8-4918-b3af-a3f81bbe519c", dataset.UniqueId());
    EXPECT_EQ("2.3.0",                                dataset.Version());
    EXPECT_EQ("http://pacificbiosciences.com/PacBioDataModel.xsd", dataset.Attribute("xmlns"));
    EXPECT_EQ("http://www.w3.org/2001/XMLSchema-instance",         dataset.Attribute("xmlns:xsi"));

    const ExternalResources& resources = dataset.ExternalResources();
    EXPECT_EQ(2, resources.Size());

    const ExternalResource& resource1 = resources[0];
    EXPECT_EQ("Third Alignments BAM",                      resource1.Name());
    EXPECT_EQ("Points to an example Alignments BAM file.", resource1.Description());
    EXPECT_EQ("AlignmentFile.AlignmentBamFile",            resource1.MetaType());
    EXPECT_EQ("file:/mnt/path/to/alignments2.bam",         resource1.ResourceId());
    EXPECT_EQ("Example",                                   resource1.Tags());
    const FileIndices& fileIndices1 = resource1.FileIndices();
    EXPECT_EQ(1, fileIndices1.Size());
    const FileIndex& pbi1 = fileIndices1[0];
    EXPECT_EQ("PacBio.Index.PacBioIndex",          pbi1.MetaType());
    EXPECT_EQ("file:/mnt/path/to/alignments2.pbi", pbi1.ResourceId());

    const ExternalResource& resource2 = resources[1];
    EXPECT_EQ("Fourth Alignments BAM",                     resource2.Name());
    EXPECT_EQ("Points to another example Alignments BAM file, by relative path.", resource2.Description());
    EXPECT_EQ("AlignmentFile.AlignmentBamFile",            resource2.MetaType());
    EXPECT_EQ("file:./alignments3.bam",                    resource2.ResourceId());
    EXPECT_EQ("Example",                                   resource2.Tags());
    const FileIndices& fileIndices2 = resource2.FileIndices();
    EXPECT_EQ(1, fileIndices2.Size());
    const FileIndex& pbi2 = fileIndices2[0];
    EXPECT_EQ("PacBio.Index.PacBioIndex",          pbi2.MetaType());
    EXPECT_EQ("file:/mnt/path/to/alignments3.pbi", pbi2.ResourceId());

    const SubDataSets& subDatasets = dataset.SubDataSets();
    EXPECT_EQ(2, subDatasets.Size());

    const DataSetBase& sub1 = subDatasets[0];
    EXPECT_EQ("HighQuality Read Alignments",          sub1.Name());
    EXPECT_EQ("ab95d0a3-94b8-4918-b3af-a3f81bbe519c", sub1.UniqueId());
    EXPECT_EQ("2.3.0",                                sub1.Version());
    const Filters& sub1Filters = sub1.Filters();
    EXPECT_EQ(1, sub1Filters.Size());
    const Filter& sub1Filter = sub1Filters[0];
    EXPECT_EQ(1, sub1Filter.Properties().Size());
    const Property& property1 = sub1Filter.Properties()[0];
    EXPECT_EQ("rq",   property1.Name());
    EXPECT_EQ(">",    property1.Operator());
    EXPECT_EQ("0.85", property1.Value());

    const DataSetBase& sub2 = subDatasets[1];
    EXPECT_EQ("Alignments to chromosome 1",          sub2.Name());
    EXPECT_EQ("ac95d0a3-94b8-4918-b3af-a3f81bbe519c", sub2.UniqueId());
    EXPECT_EQ("2.3.0",                                sub2.Version());
    const Filters& sub2Filters = sub2.Filters();
    EXPECT_EQ(1, sub2Filters.Size());
    const Filter& sub2Filter = sub2Filters[0];
    EXPECT_EQ(1, sub2Filter.Properties().Size());
    const Property& property2 = sub2Filter.Properties()[0];
    EXPECT_EQ("RNAME",   property2.Name());
    EXPECT_EQ("==",    property2.Operator());
    EXPECT_EQ("chr1", property2.Value());
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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
    EXPECT_EQ(string("http://www.w3.org/2001/XMLSchema-instance"),        dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xsi:schemaLocation"));

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
    EXPECT_EQ(string("http://pacificbiosciences.com/PacBioDatasets.xsd"), dataset.Attribute("xmlns"));
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

TEST(DataSetIOTest, InspectMalformedXml)
{
    const string xmlFn = tests::Data_Dir + "/dataset/malformed.xml";

    DataSet ds(xmlFn);
    stringstream s;
    ds.SaveToStream(s);

    const string expected =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<SubreadSet Description=\"Merged dataset from 1 files using DatasetMerger 0.1.2\" "
                    "MetaType=\"PacBio.DataSet.HdfSubreadSet\" Name=\"Subreads from runr000013_42267_150403\" "
                    "Tags=\"pacbio.secondary.instrument=RS\" TimeStampedName=\"hdfsubreadset_2015-08-19T15:39:36.331-07:00\" "
                    "UniqueId=\"b4741521-2a4c-42df-8a13-0a755ca9ed1e\" Version=\"0.5\" "
                    "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
                    "xmlns:ns0=\"http://pacificbiosciences.com/PacBioBaseDataModel.xsd\" "
                    "xmlns:ns1=\"http://pacificbiosciences.com/PacBioSampleInfo.xsd\" "
                    "xmlns:ns2=\"http://pacificbiosciences.com/PacBioCollectionMetadata.xsd\" "
                    "xmlns:ns3=\"http://pacificbiosciences.com/PacBioReagentKit.xsd\" "
                    "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                    "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t<ns0:ExternalResources>\n"
        "\t\t<ns0:ExternalResource MetaType=\"SubreadFile.SubreadBamFile\" "
                                  "ResourceId=\"file:///mnt/secondary-siv/jenkins/jenkins-bot01/workspace/Ubuntu1404_Mainline_SA3_Tiny_tests/software/smrtanalysis/siv/testkit-jobs/sa3_pipelines/mapping/tiny/job_output-ubuntu1404/tasks/pbsmrtpipe.tasks.h5_subreads_to_subread-0//mnt/secondary-siv/jenkins/jenkins-bot01/workspace/Ubuntu1404_Mainline_SA3_Tiny_tests/software/smrtanalysis/siv/testkit-jobs/sa3_pipelines/mapping/tiny/job_output-ubuntu1404/tasks/pbsmrtpipe.tasks.h5_subreads_to_subread-0/file.subreads.subreads.bam\" "
                                  "TimeStampedName=\"SubreadFile.SubreadBamFile_00000000000000\" "
                                  "UniqueId=\"251acf71-9eb0-489e-9dd1-cdbd11432753\" />\n"
        "\t</ns0:ExternalResources>\n"
        "\t<DataSetMetadata>\n"
        "\t\t<TotalLength>50000000</TotalLength>\n"
        "\t\t<NumRecords>150000</NumRecords>\n"
        "\t\t<ns2:Collections>\n"
        "\t\t\t<ns2:CollectionMetadata Context=\"m150404_101626_42267_c100807920800000001823174110291514_s1_p0\" "
                                      "InstrumentId=\"1\" InstrumentName=\"42267\" MetaType=\"PacBio.Collection\" "
                                      "TimeStampedName=\"m150404_101626_42267_c100807920800000001823174110291514_s1_p0\" "
                                      "UniqueId=\"d66c8372-2b70-4dcf-b64f-9f8b5cc351fd\">\n"
        "\t\t\t\t<ns2:InstCtrlVer>2.3.0.1.142990</ns2:InstCtrlVer>\n"
        "\t\t\t\t<ns2:SigProcVer>NRT@172.31.128.10:8082, SwVer=2301.142990, HwVer=1.0</ns2:SigProcVer>\n"
        "\t\t\t\t<ns2:RunDetails>\n"
        "\t\t\t\t\t<ns2:RunId>r000013_42267_150403</ns2:RunId>\n"
        "\t\t\t\t\t<ns2:Name>Inst42267-040315-SAT-100pM-2kb-P6C4</ns2:Name>\n"
        "\t\t\t\t</ns2:RunDetails>\n"
        "\t\t\t\t<ns2:WellSample Name=\"Inst42267-040315-SAT-100pM-2kb-P6C4\">\n"
        "\t\t\t\t\t<ns2:PlateId>Inst42267-040315-SAT-100pM-2kb-P6C4</ns2:PlateId>\n"
        "\t\t\t\t\t<ns2:WellName>Inst42267-040315-SAT-100pM-2kb-P6C4</ns2:WellName>\n"
        "\t\t\t\t\t<ns2:Concentration>0.0</ns2:Concentration>\n"
        "\t\t\t\t\t<ns2:SampleReuseEnabled>false</ns2:SampleReuseEnabled>\n"
        "\t\t\t\t\t<ns2:StageHotstartEnabled>false</ns2:StageHotstartEnabled>\n"
        "\t\t\t\t\t<ns2:SizeSelectionEnabled>false</ns2:SizeSelectionEnabled>\n"
        "\t\t\t\t\t<ns2:UseCount>1</ns2:UseCount>\n"
        "\t\t\t\t\t<ns1:BioSamplePointers>\n"
        "\t\t\t\t\t\t<ns1:BioSamplePointer>251acf71-9eb0-489e-9dd1-cdbd11432752</ns1:BioSamplePointer>\n"
        "\t\t\t\t\t</ns1:BioSamplePointers>\n"
        "\t\t\t\t</ns2:WellSample>\n"
        "\t\t\t\t<ns2:Automation>\n"
        "\t\t\t\t\t<ns0:AutomationParameters>\n"
        "\t\t\t\t\t\t<ns0:AutomationParameter />\n"
        "\t\t\t\t\t</ns0:AutomationParameters>\n"
        "\t\t\t\t</ns2:Automation>\n"
        "\t\t\t\t<ns2:CollectionNumber>7</ns2:CollectionNumber>\n"
        "\t\t\t\t<ns2:CellIndex>4</ns2:CellIndex>\n"
        "\t\t\t\t<ns2:CellPac Barcode=\"10080792080000000182317411029151\" />\n"
        "\t\t\t\t<ns2:Primary>\n"
        "\t\t\t\t\t<ns2:AutomationName>BasecallerV1</ns2:AutomationName>\n"
        "\t\t\t\t\t<ns2:ConfigFileName>2-3-0_P6-C4.xml</ns2:ConfigFileName>\n"
        "\t\t\t\t\t<ns2:SequencingCondition />\n"
        "\t\t\t\t\t<ns2:OutputOptions>\n"
        "\t\t\t\t\t\t<ns2:ResultsFolder>Analysis_Results</ns2:ResultsFolder>\n"
        "\t\t\t\t\t\t<ns2:CollectionPathUri>rsy://mp-rsync/vol55//RS_DATA_STAGING/42267/Inst42267-040315-SAT-100pM-2kb-P6C4_13/A04_7/</ns2:CollectionPathUri>\n"
        "\t\t\t\t\t\t<ns2:CopyFiles>\n"
        "\t\t\t\t\t\t\t<ns2:CollectionFileCopy>Fasta</ns2:CollectionFileCopy>\n"
        "\t\t\t\t\t\t</ns2:CopyFiles>\n"
        "\t\t\t\t\t\t<ns2:Readout>Bases</ns2:Readout>\n"
        "\t\t\t\t\t\t<ns2:MetricsVerbosity>Minimal</ns2:MetricsVerbosity>\n"
        "\t\t\t\t\t</ns2:OutputOptions>\n"
        "\t\t\t\t</ns2:Primary>\n"
        "\t\t\t</ns2:CollectionMetadata>\n"
        "\t\t</ns2:Collections>\n"
        "\t\t<ns1:BioSamples>\n"
        "\t\t\t<ns1:BioSample Description=\"Inst42267-SAT-100pM-2kbLambda-P6C4-Std120_CPS_040315\" "
                            "MetaType=\"PacBio.Sample\" Name=\"Inst42267-040315-SAT-100pM-2kb-P6C4\" "
                            "TimeStampedName=\"biosample_2015-08-19T15:39:36.331-07:00\" UniqueId=\"251acf71-9eb0-489e-9dd1-cdbd11432752\" />\n"
        "\t\t</ns1:BioSamples>\n"
        "\t</DataSetMetadata>\n"
        "</SubreadSet>\n";

    EXPECT_EQ(expected, s.str());
}

TEST(DataSetIOTest, RelativePathCarriedThroughOk_FromString)
{
    const string inputXml =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:AlignmentSet "
            "CreatedAt=\"2015-01-27T09:00:01\" "
            "MetaType=\"PacBio.DataSet.AlignmentSet\" "
            "Name=\"DataSet_AlignmentSet\" "
            "Tags=\"barcode moreTags mapping mytags\" "
            "TimeStampedName=\"biosample_2015-08-19T15:39:36.331-07:00\" "
            "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" "
            "Version=\"2.3.0\" "
            "xmlns=\"http://pacificbiosciences.com/PacBioDataModel.xsd\" "
            "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
            "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
            "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDataModel.xsd\">\n"
        "\t<pbbase:ExternalResources>\n"
        "\t\t<pbbase:ExternalResource "
                "Description=\"Points to an example Alignments BAM file.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" "
                "Name=\"Third Alignments BAM\" "
                "ResourceId=\"../path/to/resource1.bam\" "
                "Tags=\"Example\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
                    "MetaType=\"PacBio.Index.PacBioIndex\" "
                    "ResourceId=\"../path/to/resource1.bam.pbi\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t\t<pbbase:ExternalResource "
                "Description=\"Points to another example Alignments BAM file, by relative path.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" "
                "Name=\"Fourth Alignments BAM\" "
                "ResourceId=\"../path/to/resource2.bam\" "
                "Tags=\"Example\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
                    "MetaType=\"PacBio.Index.PacBioIndex\" "
                    "ResourceId=\"../path/to/resource2.bam.pbi\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t</pbbase:ExternalResources>\n"
        "</pbds:AlignmentSet>\n";

    auto dataset = DataSet::FromXml(inputXml);

    stringstream stream;
    dataset.SaveToStream(stream);
    auto outputXml = stream.str();

    EXPECT_EQ(inputXml, outputXml);
}

TEST(DataSetIOTest, RelativePathCarriedThroughOk_FromFile)
{
    DataSet dataset(tests::Data_Dir + "/relative/relative.xml");
    auto resources = dataset.ExternalResources();
    EXPECT_EQ("./a/test.bam",  resources[0].ResourceId());
    EXPECT_EQ("./b/test1.bam", resources[1].ResourceId());
    EXPECT_EQ("./b/test2.bam", resources[2].ResourceId());

    stringstream out;
    dataset.SaveToStream(out);

    auto newDataset = DataSet::FromXml(out.str());
    auto newResources = newDataset.ExternalResources();
    EXPECT_EQ("./a/test.bam",  newResources[0].ResourceId());
    EXPECT_EQ("./b/test1.bam", newResources[1].ResourceId());
    EXPECT_EQ("./b/test2.bam", newResources[2].ResourceId());
}

TEST(DataSetIOTest, DataSetFromRelativeBamFilename)
{
    // cache initial directory and move to location so we can test relatvie filename ok
    const string startingDirectory = internal::FileUtils::CurrentWorkingDirectory();

    const string targetDirectory = tests::Data_Dir + "/dataset";
    changeCurrentDirectory(targetDirectory);
    ASSERT_EQ(targetDirectory, internal::FileUtils::CurrentWorkingDirectory());

    EXPECT_NO_THROW(
    {
        const string relativeBamFn = "../phi29.bam";
        const DataSet ds(relativeBamFn);
        const auto& files = ds.BamFiles();
        EXPECT_EQ(1, files.size());
    });

    // restore working directory
    changeCurrentDirectory(startingDirectory);
}


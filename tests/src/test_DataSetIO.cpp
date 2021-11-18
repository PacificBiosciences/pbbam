#include <pbbam/DataSet.h>

#include <cstddef>

#include <fstream>

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include <gtest/gtest.h>
#include <boost/algorithm/string.hpp>

#include <pbbam/internal/DataSetElement.h>

#include "../src/FileUtils.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace DataSetIOTests {

const std::string alignedBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const std::string bamGroupFofn = PbbamTestsConfig::Generated_Dir + "/group.fofn";

const std::string ali1XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/ali1.xml";
const std::string ali2XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/ali2.xml";
const std::string ali3XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/ali3.xml";
const std::string ali4XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/ali4.xml";
const std::string mappingStaggeredXmlFn =
    PbbamTestsConfig::Data_Dir + "/dataset/bam_mapping_staggered.xml";
const std::string barcodeXmlFn = PbbamTestsConfig::Data_Dir + "/dataset/barcode.dataset.xml";
const std::string ccsReadXmlFn = PbbamTestsConfig::Data_Dir + "/dataset/ccsread.dataset.xml";
const std::string lambdaContigsXmlFn = PbbamTestsConfig::Data_Dir + "/dataset/lambda_contigs.xml";
const std::string pbalchemyXmlFn = PbbamTestsConfig::Data_Dir + "/dataset/pbalchemy10kbp.xml";
const std::string referenceXmlFn = PbbamTestsConfig::Data_Dir + "/dataset/reference.dataset.xml";
const std::string subread1XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/subread_dataset1.xml";
const std::string subread2XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/subread_dataset2.xml";
const std::string subread3XmlFn = PbbamTestsConfig::Data_Dir + "/dataset/subread_dataset3.xml";
const std::string transformedXmlFn =
    PbbamTestsConfig::Data_Dir + "/dataset/transformed_rs_subread_dataset.xml";

static void TestFromXmlString();
static void TestAli1Xml();
static void TestAli2Xml();
static void TestAli3Xml();
static void TestAli4Xml();
static void TestMappingStaggeredXml();
static void TestBarcodeXml();
static void TestCcsReadXml();
static void TestLambdaContigsXml();
static void TestPbalchemyXml();
static void TestReferenceXml();
static void TestSubread1Xml();
static void TestSubread2Xml();
static void TestSubread3Xml();
static void TestTransformedXml();

static void changeCurrentDirectory(const std::string& dir) { ASSERT_EQ(0, chdir(dir.c_str())); }

}  // namespace DataSetIOTests

TEST(BAM_DataSetIO, can_create_from_single_bam_path)
{
    DataSet dataset(DataSetIOTests::alignedBamFn);

    EXPECT_EQ(1, dataset.ExternalResources().Size());
    const ExternalResource& bamRef = dataset.ExternalResources()[0];

    EXPECT_EQ(DataSetIOTests::alignedBamFn, bamRef.ResourceId());
}

TEST(BAM_DataSetIO, can_create_from_multiple_bam_paths)
{
    std::ifstream fofn(DataSetIOTests::bamGroupFofn);
    std::vector<std::string> files;
    std::string file;
    while (std::getline(fofn, file)) {
        if (!file.empty()) {
            files.emplace_back(file);
        }
    }
    DataSet dataset(files);
    EXPECT_EQ(3, dataset.ExternalResources().Size());
}

TEST(BAM_DataSetIO, can_create_from_bam_file_object)
{
    BamFile bamFile(DataSetIOTests::alignedBamFn);
    DataSet dataset(bamFile.Filename());

    EXPECT_EQ(1, dataset.ExternalResources().Size());
    const ExternalResource& bamRef = dataset.ExternalResources()[0];

    EXPECT_EQ(DataSetIOTests::alignedBamFn, bamRef.ResourceId());
}

TEST(BAM_DataSetIO, FromFofn)
{
    DataSet dataset(DataSetIOTests::bamGroupFofn);
    EXPECT_EQ(3, dataset.ExternalResources().Size());
}

TEST(BAM_DataSetIO, can_create_from_fofn) { EXPECT_NO_THROW(DataSetIOTests::TestFromXmlString()); }

TEST(BAM_DataSetIO, can_create_from_xml)
{
    EXPECT_NO_THROW(DataSetIOTests::TestAli1Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestAli2Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestAli3Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestAli4Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestMappingStaggeredXml());
    EXPECT_NO_THROW(DataSetIOTests::TestBarcodeXml());
    EXPECT_NO_THROW(DataSetIOTests::TestCcsReadXml());
    EXPECT_NO_THROW(DataSetIOTests::TestLambdaContigsXml());
    EXPECT_NO_THROW(DataSetIOTests::TestPbalchemyXml());
    EXPECT_NO_THROW(DataSetIOTests::TestReferenceXml());
    EXPECT_NO_THROW(DataSetIOTests::TestSubread1Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestSubread2Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestSubread3Xml());
    EXPECT_NO_THROW(DataSetIOTests::TestTransformedXml());
}

TEST(BAM_DataSetIO, throws_on_nonexistent_fofn)
{
    EXPECT_THROW(DataSet("does/not/exist.fofn"), std::exception);
}

TEST(BAM_DataSetIO, throws_on_nonexistent_xml)
{
    EXPECT_THROW(DataSet("does/not/exist.xml"), std::exception);
}

TEST(BAM_DataSetIO, throws_on_supported_extension)
{
    EXPECT_THROW(DataSet("bad/extension.foo"), std::exception);
}

TEST(BAM_DataSetIO, throws_if_cannot_save_to_file)
{
    DataSet ds;
    EXPECT_THROW(ds.Save("fake_directory_that_should_not_exist/out.xml"), std::exception);
}

TEST(BAM_DataSetIO, can_write_normal_alignmentset_as_xml)
{
    // top-level data
    DataSet dataset(DataSet::ALIGNMENT);
    dataset.CreatedAt("2015-01-27T09:00:01");
    dataset.MetaType("PacBio.DataSet.AlignmentSet");
    dataset.Name("DataSet_AlignmentSet");
    dataset.Tags("barcode moreTags mapping mytags");
    dataset.TimeStampedName("my_tsn");
    dataset.UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c");
    dataset.Attribute("xmlns", "http://pacificbiosciences.com/PacBioDatasets.xsd")
        .Attribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
        .Attribute("xsi:schemaLocation", "http://pacificbiosciences.com/PacBioDatasets.xsd");

    // external resources
    ExternalResource resource1("AlignmentFile.AlignmentBamFile", "/mnt/path/to/alignments2.bam");
    resource1.Name("Third Alignments BAM");
    resource1.CreatedAt("2015-01-27T09:00:01");
    resource1.Description("Points to an example Alignments BAM file.");
    resource1.Tags("Example");
    resource1.TimeStampedName("my_tsn");
    resource1.UniqueId("my_uuid");
    FileIndex pbi1("PacBio.Index.PacBioIndex", "/mnt/path/to/alignments2.pbi");
    pbi1.CreatedAt("2015-01-27T09:00:01");
    pbi1.TimeStampedName("my_tsn");
    pbi1.UniqueId("my_uuid");
    resource1.FileIndices().Add(pbi1);
    dataset.ExternalResources().Add(resource1);

    ExternalResource resource2("AlignmentFile.AlignmentBamFile", "./alignments3.bam");
    resource2.CreatedAt("2015-01-27T09:00:01");
    resource2.Name("Fourth Alignments BAM");
    resource2.Description("Points to another example Alignments BAM file, by relative path.");
    resource2.Tags("Example");
    resource2.TimeStampedName("my_tsn");
    resource2.UniqueId("my_uuid");
    FileIndex pbi2("PacBio.Index.PacBioIndex", "./alignments3.pbi");
    pbi2.CreatedAt("2015-01-27T09:00:01");
    pbi2.TimeStampedName("my_tsn");
    pbi2.UniqueId("my_uuid");

    resource2.FileIndices().Add(pbi2);
    dataset.ExternalResources().Add(resource2);

    // sub-datasets with filters
    DataSetBase subDataSet1;
    subDataSet1.CreatedAt("2015-01-27T09:00:01");
    subDataSet1.Name("HighQuality Read Alignments");
    subDataSet1.TimeStampedName("my_tsn");
    subDataSet1.UniqueId("ab95d0a3-94b8-4918-b3af-a3f81bbe519c");
    Filter filter1;
    filter1.Properties().Add(Property("rq", "0.85", ">"));
    subDataSet1.Filters().Add(filter1);
    dataset.SubDataSets().Add(subDataSet1);

    DataSetBase subDataSet2;
    subDataSet2.CreatedAt("2015-01-27T09:00:01");
    subDataSet2.Name("Alignments to chromosome 1");
    subDataSet2.TimeStampedName("my_tsn");
    subDataSet2.UniqueId("ac95d0a3-94b8-4918-b3af-a3f81bbe519c");
    Filter filter2;
    filter2.Properties().Add(Property("RNAME", "chr1", "=="));
    subDataSet2.Filters().Add(filter2);
    dataset.SubDataSets().Add(subDataSet2);

    // write dataset
    const std::string expectedXml{
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
        "CreatedAt=\"2015-01-27T09:00:01\" "
        "Description=\"Points to an example Alignments BAM file.\" "
        "MetaType=\"AlignmentFile.AlignmentBamFile\" "
        "Name=\"Third Alignments BAM\" "
        "ResourceId=\"/mnt/path/to/alignments2.bam\" "
        "Tags=\"Example\" "
        "TimeStampedName=\"my_tsn\" "
        "UniqueId=\"my_uuid\" Version=\"3.0.1\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
        "CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.Index.PacBioIndex\" "
        "ResourceId=\"/mnt/path/to/alignments2.pbi\" "
        "TimeStampedName=\"my_tsn\" "
        "UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t\t<pbbase:ExternalResource "
        "CreatedAt=\"2015-01-27T09:00:01\" "
        "Description=\"Points to another example Alignments BAM file, by relative path.\" "
        "MetaType=\"AlignmentFile.AlignmentBamFile\" "
        "Name=\"Fourth Alignments BAM\" "
        "ResourceId=\"" +
        dataset.Path() +
        "/alignments3.bam\" "
        "Tags=\"Example\" "
        "TimeStampedName=\"my_tsn\" "
        "UniqueId=\"my_uuid\" Version=\"3.0.1\">\n"
        "\t\t\t<pbbase:FileIndices>\n"
        "\t\t\t\t<pbbase:FileIndex "
        "CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.Index.PacBioIndex\" "
        "ResourceId=\"" +
        dataset.Path() +
        "/alignments3.pbi\" "
        "TimeStampedName=\"my_tsn\" "
        "UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t\t\t</pbbase:FileIndices>\n"
        "\t\t</pbbase:ExternalResource>\n"
        "\t</pbbase:ExternalResources>\n"
        "\t<pbds:DataSets>\n"
        "\t\t<pbds:DataSet "
        "CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.DataSet.DataSet\" "
        "Name=\"HighQuality Read Alignments\" "
        "TimeStampedName=\"my_tsn\" "
        "UniqueId=\"ab95d0a3-94b8-4918-b3af-a3f81bbe519c\" "
        "Version=\"3.0.1\">\n"
        "\t\t\t<pbds:Filters>\n"
        "\t\t\t\t<pbds:Filter>\n"
        "\t\t\t\t\t<pbbase:Properties>\n"
        "\t\t\t\t\t\t<pbbase:Property Name=\"rq\" Operator=\"&gt;\" Value=\"0.85\" />\n"
        "\t\t\t\t\t</pbbase:Properties>\n"
        "\t\t\t\t</pbds:Filter>\n"
        "\t\t\t</pbds:Filters>\n"
        "\t\t</pbds:DataSet>\n"
        "\t\t<pbds:DataSet "
        "CreatedAt=\"2015-01-27T09:00:01\" "
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
        "</pbds:AlignmentSet>\n"};

    std::ostringstream s;
    dataset.SaveToStream(s);
    EXPECT_EQ(expectedXml, s.str());
}

TEST(BAM_DataSetIO, can_write_normal_contigset_as_xml)
{
    // top-level data
    ContigSet dataset;
    dataset.CreatedAt("2015-01-27T09:00:01");
    dataset.Name("DataSet_ContigSet");
    dataset.Tags("barcode moreTags mapping mytags");
    dataset.TimeStampedName("my_tsn");
    dataset.UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c");

    // write dataset
    const std::string expectedXml{
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:ContigSet "
        "CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.DataSet.ContigSet\" "
        "Name=\"DataSet_ContigSet\" "
        "Tags=\"barcode moreTags mapping mytags\" "
        "TimeStampedName=\"my_tsn\" "
        "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
        "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" />\n"};

    std::ostringstream s;
    dataset.SaveToStream(s);
    EXPECT_EQ(expectedXml, s.str());
}

namespace DataSetIOTests {

static void TestFromXmlString()
{
    const std::string inputXml{
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
        "</pbds:AlignmentSet>\n"};

    const DataSet dataset = DataSet::FromXml(inputXml);

    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ("2015-01-27T09:00:01", dataset.CreatedAt());
    EXPECT_EQ("PacBio.DataSet.AlignmentSet", dataset.MetaType());
    EXPECT_EQ("DataSet_AlignmentSet", dataset.Name());
    EXPECT_EQ("barcode moreTags mapping mytags", dataset.Tags());
    EXPECT_EQ("b095d0a3-94b8-4918-b3af-a3f81bbe519c", dataset.UniqueId());
    EXPECT_EQ("2.3.0", dataset.Version());
    EXPECT_EQ("http://pacificbiosciences.com/PacBioDataModel.xsd", dataset.Attribute("xmlns"));
    EXPECT_EQ("http://www.w3.org/2001/XMLSchema-instance", dataset.Attribute("xmlns:xsi"));

    const ExternalResources& resources = dataset.ExternalResources();
    EXPECT_EQ(2, resources.NumChildren());

    const ExternalResource& resource1 = resources[0];
    EXPECT_EQ("Third Alignments BAM", resource1.Name());
    EXPECT_EQ("Points to an example Alignments BAM file.", resource1.Description());
    EXPECT_EQ("AlignmentFile.AlignmentBamFile", resource1.MetaType());
    EXPECT_EQ("file:/mnt/path/to/alignments2.bam", resource1.ResourceId());
    EXPECT_EQ("Example", resource1.Tags());
    const FileIndices& fileIndices1 = resource1.FileIndices();
    EXPECT_EQ(1, fileIndices1.Size());
    const FileIndex& pbi1 = fileIndices1[0];
    EXPECT_EQ("PacBio.Index.PacBioIndex", pbi1.MetaType());
    EXPECT_EQ("file:/mnt/path/to/alignments2.pbi", pbi1.ResourceId());

    const ExternalResource& resource2 = resources[1];
    EXPECT_EQ("Fourth Alignments BAM", resource2.Name());
    EXPECT_EQ("Points to another example Alignments BAM file, by relative path.",
              resource2.Description());
    EXPECT_EQ("AlignmentFile.AlignmentBamFile", resource2.MetaType());
    EXPECT_EQ("file:./alignments3.bam", resource2.ResourceId());
    EXPECT_EQ("Example", resource2.Tags());
    const FileIndices& fileIndices2 = resource2.FileIndices();
    EXPECT_EQ(1, fileIndices2.Size());
    const FileIndex& pbi2 = fileIndices2[0];
    EXPECT_EQ("PacBio.Index.PacBioIndex", pbi2.MetaType());
    EXPECT_EQ("file:/mnt/path/to/alignments3.pbi", pbi2.ResourceId());

    const SubDataSets& subDatasets = dataset.SubDataSets();
    EXPECT_EQ(2, subDatasets.Size());

    const DataSetBase& sub1 = subDatasets[0];
    EXPECT_EQ("HighQuality Read Alignments", sub1.Name());
    EXPECT_EQ("ab95d0a3-94b8-4918-b3af-a3f81bbe519c", sub1.UniqueId());
    EXPECT_EQ("2.3.0", sub1.Version());
    const Filters& sub1Filters = sub1.Filters();
    EXPECT_EQ(1, sub1Filters.Size());
    const Filter& sub1Filter = sub1Filters[0];
    EXPECT_EQ(1, sub1Filter.Properties().Size());
    const Property& property1 = sub1Filter.Properties()[0];
    EXPECT_EQ("rq", property1.Name());
    EXPECT_EQ(">", property1.Operator());
    EXPECT_EQ("0.85", property1.Value());

    const DataSetBase& sub2 = subDatasets[1];
    EXPECT_EQ("Alignments to chromosome 1", sub2.Name());
    EXPECT_EQ("ac95d0a3-94b8-4918-b3af-a3f81bbe519c", sub2.UniqueId());
    EXPECT_EQ("2.3.0", sub2.Version());
    const Filters& sub2Filters = sub2.Filters();
    EXPECT_EQ(1, sub2Filters.Size());
    const Filter& sub2Filter = sub2Filters[0];
    EXPECT_EQ(1, sub2Filter.Properties().Size());
    const Property& property2 = sub2Filter.Properties()[0];
    EXPECT_EQ("RNAME", property2.Name());
    EXPECT_EQ("==", property2.Operator());
    EXPECT_EQ("chr1", property2.Value());
}

static void TestAli1Xml()
{
    const DataSet dataset(DataSetIOTests::ali1XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.AlignmentSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_AlignmentSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Alignments BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments0.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments0.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(
                std::string("Points to another example Alignments BAM file, by relative path."),
                resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:./alignments1.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments1.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.85"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("RNAME"), property.Name());
            EXPECT_EQ(std::string("chr1"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }
}

static void TestAli2Xml()
{
    const DataSet dataset(DataSetIOTests::ali2XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.AlignmentSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_AlignmentSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Alignments BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments2.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments2.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(
                std::string("Points to another example Alignments BAM file, by relative path."),
                resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:./alignments3.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments3.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.85"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("RNAME"), property.Name());
            EXPECT_EQ(std::string("chr1"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }
}

static void TestAli3Xml()
{
    const DataSet dataset(DataSetIOTests::ali3XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.AlignmentSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_AlignmentSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Alignments BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments2.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments2.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(
                std::string("Points to another example Alignments BAM file, by relative path."),
                resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:./alignments3.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments3.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.75"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("RNAME"), property.Name());
            EXPECT_EQ(std::string("chr1"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }
}

static void TestAli4Xml()
{
    const DataSet dataset(DataSetIOTests::ali4XmlFn);
    EXPECT_EQ(DataSet::ALIGNMENT, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.AlignmentSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_AlignmentSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Alignments BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Alignments BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments0.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments0.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Alignments BAM"), resource.Name());
            EXPECT_EQ(
                std::string("Points to another example Alignments BAM file, by relative path."),
                resource.Description());
            EXPECT_EQ(std::string("AlignmentFile.AlignmentBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:./alignments1.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/alignments1.pbi"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("HighQuality Read Alignments"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.85"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            EXPECT_EQ(std::string(""), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string("Alignments to chromosome 1"), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const Filters& filters = subdataset.Filters();
            ASSERT_EQ(1, filters.Size());
            const Filter& filter = filters[0];
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("RNAME"), property.Name());
            EXPECT_EQ(std::string("chr1"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }
}

static void TestMappingStaggeredXml()
{
    const DataSet dataset(DataSetIOTests::mappingStaggeredXmlFn);
    EXPECT_EQ(DataSet::GENERIC, dataset.Type());
    EXPECT_EQ(std::string("2015-05-13T10:58:26"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.DataSet"), dataset.MetaType());
    EXPECT_EQ(std::string(""), dataset.Name());
    EXPECT_EQ(std::string(""), dataset.Tags());
    EXPECT_EQ(std::string("30f72098-bc5b-e06b-566c-8b28dda909a8"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string(""), resource.Name());
            EXPECT_EQ(std::string(""), resource.Description());
            EXPECT_EQ(std::string(""), resource.MetaType());
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_1.bam"), resource.ResourceId());
            EXPECT_EQ(std::string(""), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_1.bam.bai"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string(""), resource.Name());
            EXPECT_EQ(std::string(""), resource.Description());
            EXPECT_EQ(std::string(""), resource.MetaType());
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_2.bam"), resource.ResourceId());
            EXPECT_EQ(std::string(""), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_2.bam.bai"), index.ResourceId());
        }
    }

    const SubDataSets& subdatasets = dataset.SubDataSets();
    ASSERT_EQ(2, subdatasets.Size());
    for (size_t i = 0; i < subdatasets.Size(); ++i) {
        const DataSetBase& subdataset = subdatasets[i];
        if (i == 0) {
            EXPECT_EQ(std::string("2015-05-13T10:58:26"), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string(""), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("c5402d06-4643-057c-e300-fe229b4e8909"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const ExternalResources& subResources = subdataset.ExternalResources();
            ASSERT_EQ(1, subResources.Size());
            const ExternalResource& resource = subResources[0];
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_2.bam"), resource.ResourceId());
            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_2.bam.bai"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("2015-05-13T10:58:26"), subdataset.CreatedAt());
            EXPECT_EQ(std::string(""), subdataset.MetaType());
            EXPECT_EQ(std::string(""), subdataset.Name());
            EXPECT_EQ(std::string(""), subdataset.Tags());
            EXPECT_EQ(std::string("f8b54a55-5fb7-706f-ab35-39afc9c86924"), subdataset.UniqueId());
            EXPECT_EQ(std::string("2.3.0"), subdataset.Version());

            const ExternalResources& subResources = subdataset.ExternalResources();
            ASSERT_EQ(1, subResources.Size());
            const ExternalResource& resource = subResources[0];
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_1.bam"), resource.ResourceId());
            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:tests/data/bam_mapping_1.bam.bai"), index.ResourceId());
        }
    }
}

static void TestBarcodeXml()
{
    const DataSet dataset(DataSetIOTests::barcodeXmlFn);
    EXPECT_EQ(DataSet::BARCODE, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.BarcodeSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_BarcodeSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(std::string("First Barcodes FASTA"), resource.Name());
    EXPECT_EQ(std::string("Points to an example Barcodes FASTA file."), resource.Description());
    EXPECT_EQ(std::string("BarcodeFile.BarcodeFastaFile"), resource.MetaType());
    EXPECT_EQ(std::string("file:///mnt/path/to/barcode.fasta"), resource.ResourceId());
    EXPECT_EQ(std::string("Example"), resource.Tags());

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(std::string("30"), metadata.NumRecords());
    EXPECT_EQ(std::string("400"), metadata.TotalLength());

    // access metadata extensions directly for now
    EXPECT_EQ(std::string("paired"), metadata.ChildText("BarcodeConstruction"));
}

static void TestCcsReadXml()
{
    const DataSet dataset(DataSetIOTests::ccsReadXmlFn);
    EXPECT_EQ(DataSet::CONSENSUS_READ, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.ConsensusReadSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_ConsensusReadSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First ConsensusRead BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example ConsensusRead BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("PacBio.ConsensusReadFile.ConsensusReadBamFile"),
                      resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/ccsreads0.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("PacBio.Index.PacBioIndex"), index.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/ccsreads0.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second ConsensusRead BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to another example ConsensusRead BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("PacBio.ConsensusReadFile.ConsensusReadBamFile"),
                      resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/ccsreads1.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("PacBio.Index.PacBioIndex"), index.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/ccsreads0.pbi"), index.ResourceId());
        }
    }
}

static void TestLambdaContigsXml()
{
    const DataSet dataset(DataSetIOTests::lambdaContigsXmlFn);
    EXPECT_EQ(DataSet::REFERENCE, dataset.Type());
    EXPECT_EQ(std::string("2015-05-28T10:56:36"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.ReferenceSet"), dataset.MetaType());
    EXPECT_EQ(std::string(""), dataset.Name());
    EXPECT_EQ(std::string(""), dataset.Tags());
    EXPECT_EQ(std::string("596e87db-34f9-d2fd-c905-b017543170e1"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(std::string("file:tests/data/lambda_contigs.fasta"), resource.ResourceId());
}

static void TestPbalchemyXml()
{
    const DataSet dataset(DataSetIOTests::pbalchemyXmlFn);
    EXPECT_EQ(DataSet::GENERIC, dataset.Type());
    EXPECT_EQ(std::string("2015-05-22T16:56:16"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.DataSet"), dataset.MetaType());
    EXPECT_EQ(std::string(""), dataset.Name());
    EXPECT_EQ(std::string(""), dataset.Tags());
    EXPECT_EQ(std::string("58e3f7c5-24c1-b58b-fbd5-37de268cc2f0"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(std::string("file:tests/data/pbalchemy10kbp.pbalign.sorted.pbver1.bam"),
              resource.ResourceId());
    const FileIndices& fileIndices = resource.FileIndices();
    ASSERT_EQ(1, fileIndices.Size());
    const FileIndex& index = fileIndices[0];
    EXPECT_EQ(std::string("file:tests/data/pbalchemy10kbp.pbalign.sorted.pbver1.bam.bai"),
              index.ResourceId());

    // TYPOs: Should be Filter Properties/Property not Parameter(s)
}

static void TestReferenceXml()
{
    const DataSet dataset(DataSetIOTests::referenceXmlFn);
    EXPECT_EQ(DataSet::REFERENCE, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.ReferenceSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_ReferenceSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(1, resources.Size());
    const ExternalResource& resource = resources[0];
    EXPECT_EQ(std::string("First References FASTA"), resource.Name());
    EXPECT_EQ(std::string("Points to an example references FASTA file."), resource.Description());
    EXPECT_EQ(std::string("PacBio.ReferenceFile.ReferenceFastaFile"), resource.MetaType());
    EXPECT_EQ(std::string("file:///mnt/path/to/reference.fasta"), resource.ResourceId());
    EXPECT_EQ(std::string("Example"), resource.Tags());
    const FileIndices& fileIndices = resource.FileIndices();
    ASSERT_EQ(2, fileIndices.Size());
    for (size_t i = 0; i < fileIndices.Size(); ++i) {
        const FileIndex& index = fileIndices[i];
        if (i == 0) {
            EXPECT_EQ(std::string("PacBio.Index.SaWriterIndex"), index.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/reference.fasta.sa"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("PacBio.Index.SamIndex"), index.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/reference.fasta.fai"), index.ResourceId());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(std::string("500"), metadata.NumRecords());
    EXPECT_EQ(std::string("5000000"), metadata.TotalLength());

    // access metadata extensions directly for now
    EXPECT_EQ(std::string("Tribble"), metadata.ChildText("Organism"));
    EXPECT_EQ(std::string("Diploid"), metadata.ChildText("Ploidy"));

    const internal::DataSetElement& contigs = metadata.Child<internal::DataSetElement>("Contigs");
    ASSERT_EQ(1, contigs.NumChildren());

    const internal::DataSetElement& contig = contigs.Child<internal::DataSetElement>(0);
    EXPECT_EQ(std::string("gi|229359445|emb|AM181176.4|"), contig.Attribute("Name"));
    EXPECT_EQ(std::string("Pseudomonas fluorescens SBW25 complete genome|quiver"),
              contig.Attribute("Description"));
    EXPECT_EQ(std::string("6722109"), contig.Attribute("Length"));
    EXPECT_EQ(std::string("f627c795efad7ce0050ed42b942d408e"), contig.Attribute("Digest"));
}

static void TestSubread1Xml()
{
    const DataSet dataset(DataSetIOTests::subread1XmlFn);
    EXPECT_EQ(DataSet::SUBREAD, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_SubreadSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Subreads BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Subreads BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads0.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads0.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Subreads BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to another example Subreads BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads1.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads0.pbi"), index.ResourceId());
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
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.75"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("QNAME"), property.Name());
            EXPECT_EQ(std::string("100/0/0_100"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(std::string("500"), metadata.NumRecords());
    EXPECT_EQ(std::string("500000"), metadata.TotalLength());
}

static void TestSubread2Xml()
{
    const DataSet dataset(DataSetIOTests::subread2XmlFn);
    EXPECT_EQ(DataSet::SUBREAD, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_SubreadSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Subreads BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Subreads BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads2.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads2.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Subreads BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to another example Subreads BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads3.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads3.pbi"), index.ResourceId());
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
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.75"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("QNAME"), property.Name());
            EXPECT_EQ(std::string("100/0/0_100"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(std::string("500"), metadata.NumRecords());
    EXPECT_EQ(std::string("500000"), metadata.TotalLength());
}

static void TestSubread3Xml()
{
    const DataSet dataset(DataSetIOTests::subread3XmlFn);
    EXPECT_EQ(DataSet::SUBREAD, dataset.Type());
    EXPECT_EQ(std::string("2015-01-27T09:00:01"), dataset.CreatedAt());
    EXPECT_EQ(std::string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(std::string("DataSet_SubreadSet"), dataset.Name());
    EXPECT_EQ(std::string("barcode moreTags mapping mytags"), dataset.Tags());
    EXPECT_EQ(std::string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(std::string("2.3.0"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema-instance"),
              dataset.Attribute("xmlns:xsi"));
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xsi:schemaLocation"));

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("First Subreads BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to an example Subreads BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads2.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads2.pbi"), index.ResourceId());
        } else {
            EXPECT_EQ(std::string("Second Subreads BAM"), resource.Name());
            EXPECT_EQ(std::string("Points to another example Subreads BAM file."),
                      resource.Description());
            EXPECT_EQ(std::string("SubreadFile.SubreadBamFile"), resource.MetaType());
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads3.bam"), resource.ResourceId());
            EXPECT_EQ(std::string("Example"), resource.Tags());

            const FileIndices& fileIndices = resource.FileIndices();
            ASSERT_EQ(1, fileIndices.Size());
            const FileIndex& index = fileIndices[0];
            EXPECT_EQ(std::string("file:///mnt/path/to/subreads3.pbi"), index.ResourceId());
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
            EXPECT_EQ(std::string("rq"), property.Name());
            EXPECT_EQ(std::string("0.85"), property.Value());
            EXPECT_EQ(std::string(">"), property.Operator());
        } else {
            const Properties& properties = filter.Properties();
            ASSERT_EQ(1, properties.Size());
            const Property& property = properties[0];
            EXPECT_EQ(std::string("QNAME"), property.Name());
            EXPECT_EQ(std::string("100/0/0_100"), property.Value());
            EXPECT_EQ(std::string("=="), property.Operator());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(std::string("500"), metadata.NumRecords());
    EXPECT_EQ(std::string("500000"), metadata.TotalLength());
}

static void TestTransformedXml()
{
    const DataSet dataset(DataSetIOTests::transformedXmlFn);
    EXPECT_EQ(DataSet::HDF_SUBREAD, dataset.Type());
    EXPECT_EQ(std::string("PacBio.DataSet.SubreadSet"), dataset.MetaType());
    EXPECT_EQ(std::string("Subreads from run r001173_42129_130607"), dataset.Name());
    EXPECT_EQ(std::string("pacbio.secondary.instrument=RS"), dataset.Tags());
    EXPECT_EQ(std::string("abbc9183-b01e-4671-8c12-19efee534647"), dataset.UniqueId());
    EXPECT_EQ(std::string("0.5"), dataset.Version());
    EXPECT_EQ(std::string("http://pacificbiosciences.com/PacBioDatasets.xsd"),
              dataset.Attribute("xmlns"));
    EXPECT_EQ(std::string("http://www.w3.org/2001/XMLSchema"), dataset.Attribute("xmlns:xs"));
    EXPECT_EQ(std::string("http://www.w3.org/2005/xpath-functions"), dataset.Attribute("xmlns:fn"));
    EXPECT_EQ(std::string("java:java.util.UUID"), dataset.Attribute("xmlns:uuid"));
    EXPECT_EQ(std::string("http://whatever"), dataset.Attribute("xmlns:bax"));

    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(3, resources.Size());
    for (size_t i = 0; i < resources.Size(); ++i) {
        const ExternalResource& resource = resources[i];
        if (i == 0) {
            EXPECT_EQ(std::string("PacBio.SubreadFile.BaxFile"), resource.MetaType());
            EXPECT_EQ(
                std::string(
                    "file:///pbi/dept/secondary/siv/testdata/LIMS/2590727/0001/Analysis_Results/"
                    "m130608_033634_42129_c100515232550000001823076608221351_s1_p0.0.bax.h5"),
                resource.ResourceId());
        } else if (i == 1) {
            EXPECT_EQ(std::string("PacBio.SubreadFile.BaxFile"), resource.MetaType());
            EXPECT_EQ(
                std::string(
                    "file:///pbi/dept/secondary/siv/testdata/LIMS/2590727/0001/Analysis_Results/"
                    "m130608_033634_42129_c100515232550000001823076608221351_s1_p0.1.bax.h5"),
                resource.ResourceId());
        } else {
            EXPECT_EQ(std::string("PacBio.SubreadFile.BaxFile"), resource.MetaType());
            EXPECT_EQ(
                std::string(
                    "file:///pbi/dept/secondary/siv/testdata/LIMS/2590727/0001/Analysis_Results/"
                    "m130608_033634_42129_c100515232550000001823076608221351_s1_p0.2.bax.h5"),
                resource.ResourceId());
        }
    }

    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(std::string("150000"), metadata.NumRecords());
    EXPECT_EQ(std::string("50000000"), metadata.TotalLength());
}

}  // namespace DataSetIOTests

TEST(BAM_DataSetIO, relative_path_is_passed_through_from_input_xml)
{
    DataSet dataset(PbbamTestsConfig::Data_Dir + "/relative/relative.xml");
    auto resources = dataset.ExternalResources();
    EXPECT_EQ("./a/test.bam", resources[0].ResourceId());
    EXPECT_EQ("./b/test1.bam", resources[1].ResourceId());
    EXPECT_EQ("./b/test2.bam", resources[2].ResourceId());

    std::ostringstream out;
    dataset.SaveToStream(out);

    auto newDataset = DataSet::FromXml(out.str());
    auto newResources = newDataset.ExternalResources();
    EXPECT_EQ("./a/test.bam", newResources[0].ResourceId());
    EXPECT_EQ("./b/test1.bam", newResources[1].ResourceId());
    EXPECT_EQ("./b/test2.bam", newResources[2].ResourceId());
}

TEST(BAM_DataSetIO, relative_path_is_passed_through_from_input_bam_path)
{
    // cache initial directory and move to location so we can test relatvie filename ok
    const std::string startingDirectory = FileUtils::CurrentWorkingDirectory();

    const std::string targetDirectory = PbbamTestsConfig::Data_Dir + "/dataset";
    DataSetIOTests::changeCurrentDirectory(targetDirectory);
    ASSERT_EQ(targetDirectory, FileUtils::CurrentWorkingDirectory());

    EXPECT_NO_THROW({
        const std::string relativeBamFn = "../phi29.bam";
        const DataSet ds(relativeBamFn);
        const auto files = ds.BamFiles();
        EXPECT_EQ(1, files.size());
    });

    // restore working directory
    DataSetIOTests::changeCurrentDirectory(startingDirectory);
}

TEST(BAM_DataSetIO, can_fetch_all_file_paths)
{
    // check  BamFiles only
    EXPECT_NO_THROW({
        const DataSet dataset(PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml");
        const auto bamFiles = dataset.BamFiles();
        EXPECT_EQ(3, bamFiles.size());
    });

    // now fetch all files (original BAMs plus PBI files)
    EXPECT_NO_THROW({
        const DataSet dataset(PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml");
        const auto allFiles = dataset.AllFiles();
        EXPECT_EQ(6, allFiles.size());
    });
}

TEST(BAM_DataSetIO, correctly_orders_metadata_default_children)
{
    DataSet dataset(DataSet::ALIGNMENT);
    dataset.CreatedAt("2015-01-27T09:00:01");
    dataset.MetaType("PacBio.DataSet.AlignmentSet");
    dataset.Name("DataSet_AlignmentSet");
    dataset.Tags("barcode moreTags mapping mytags");
    dataset.TimeStampedName("my_time_stamped_name");
    dataset.UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c");
    dataset.Attribute("xmlns", "http://pacificbiosciences.com/PacBioDatasets.xsd")
        .Attribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
        .Attribute("xsi:schemaLocation", "http://pacificbiosciences.com/PacBioDatasets.xsd");

    ExternalResource ext("Fake.MetaType", "filename");
    ext.CreatedAt("2015-01-27T09:00:01");
    ext.TimeStampedName("custom_tsn").UniqueId("my_uuid");
    dataset.ExternalResources().Add(ext);

    const auto numRecords = std::to_string(42);
    const auto totalLength = std::to_string(1000);
    DataSetMetadata metadata(numRecords, totalLength);
    dataset.Metadata(metadata);

    std::ostringstream s;
    dataset.SaveToStream(s);

    const std::string result = s.str();
    const size_t xmlnsFound = result.find("xmlns=");
    const size_t xmlnsXsiFound = result.find("xmlns:xsi=");
    const size_t xsiSchemaLocationFound = result.find("xsi:schemaLocation=");
    const size_t xmlnsPbbaseFound = result.find("xmlns:pbbase=");
    const size_t xmlnsPbdsFound = result.find("xmlns:pbds=");

    EXPECT_TRUE(xmlnsFound < xmlnsXsiFound);
    EXPECT_TRUE(xmlnsXsiFound < xsiSchemaLocationFound);
    EXPECT_TRUE(xsiSchemaLocationFound < xmlnsPbbaseFound);
    EXPECT_TRUE(xmlnsPbbaseFound < xmlnsPbdsFound);
}

TEST(BAM_DataSetIO, can_make_referenceset_from_subdataset)
{
    // ReferenceSet with ReferenceSet subdataset
    const std::string referenceSetXml{
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:ReferenceSet CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.DataSet.ReferenceSet\" "
        "Name=\"DataSet_ReferenceSet\" Tags=\"barcode moreTags mapping mytags\" "
        "TimeStampedName=\"my_time_stamped_name\" "
        "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
        "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:pbbase=\"http://pacificbiosciences.com/PacBioBaseDataModel.xsd\" "
        "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t<pbbase:ExternalResources>\n"
        "\t\t<pbbase:ExternalResource MetaType=\"Fake.MetaType\" ResourceId=\"filename\" "
        "TimeStampedName=\"custom_tsn\" UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t</pbbase:ExternalResources>\n"
        "\t<pbds:DataSets>\n"
        "\t\t<pbds:ReferenceSet> CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.DataSet.ReferenceSet\" "
        "Name=\"DataSet_ReferenceSet\" Tags=\"barcode moreTags mapping mytags\" "
        "TimeStampedName=\"my_time_stamped_name\" "
        "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
        "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:pbbase=\"http://pacificbiosciences.com/PacBioBaseDataModel.xsd\" "
        "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t\t\t<pbds:DataSetMetadata>\n"
        "\t\t\t\t<pbds:TotalLength>1000</pbds:TotalLength>\n"
        "\t\t\t\t<pbds:NumRecords>42</pbds:NumRecords>\n"
        "\t\t\t</pbds:DataSetMetadata>\n"
        "\t\t</pbds:ReferenceSet>\n"
        "\t</pbds:DataSets>\n"
        "\t<pbds:DataSetMetadata>\n"
        "\t\t<pbds:TotalLength>1000</pbds:TotalLength>\n"
        "\t\t<pbds:NumRecords>42</pbds:NumRecords>\n"
        "\t</pbds:DataSetMetadata>\n"
        "</pbds:ReferenceSet>\n"};

    EXPECT_NO_THROW(DataSet::FromXml(referenceSetXml));

    // AlignmentSet with ReferenceSet subdataset
    const std::string alignmentSetXml{
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<pbds:AlignmentSet CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.DataSet.AlignmentSet\" "
        "Name=\"DataSet_AlignmentSet\" Tags=\"barcode moreTags mapping mytags\" "
        "TimeStampedName=\"my_time_stamped_name\" "
        "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
        "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:pbbase=\"http://pacificbiosciences.com/PacBioBaseDataModel.xsd\" "
        "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t<pbbase:ExternalResources>\n"
        "\t\t<pbbase:ExternalResource MetaType=\"Fake.MetaType\" ResourceId=\"filename\" "
        "TimeStampedName=\"custom_tsn\" UniqueId=\"my_uuid\" Version=\"3.0.1\" />\n"
        "\t</pbbase:ExternalResources>\n"
        "\t<pbds:DataSets>\n"
        "\t\t<pbds:ReferenceSet> CreatedAt=\"2015-01-27T09:00:01\" "
        "MetaType=\"PacBio.DataSet.ReferenceSet\" "
        "Name=\"DataSet_ReferenceSet\" Tags=\"barcode moreTags mapping mytags\" "
        "TimeStampedName=\"my_time_stamped_name\" "
        "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"3.0.1\" "
        "xmlns=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDatasets.xsd\" "
        "xmlns:pbbase=\"http://pacificbiosciences.com/PacBioBaseDataModel.xsd\" "
        "xmlns:pbds=\"http://pacificbiosciences.com/PacBioDatasets.xsd\">\n"
        "\t\t\t<pbds:DataSetMetadata>\n"
        "\t\t\t\t<pbds:TotalLength>1000</pbds:TotalLength>\n"
        "\t\t\t\t<pbds:NumRecords>42</pbds:NumRecords>\n"
        "\t\t\t</pbds:DataSetMetadata>\n"
        "\t\t</pbds:ReferenceSet>\n"
        "\t</pbds:DataSets>\n"
        "\t<pbds:DataSetMetadata>\n"
        "\t\t<pbds:TotalLength>1000</pbds:TotalLength>\n"
        "\t\t<pbds:NumRecords>42</pbds:NumRecords>\n"
        "\t</pbds:DataSetMetadata>\n"
        "</pbds:AlignmentSet>\n"};

    EXPECT_NO_THROW(DataSet::FromXml(alignmentSetXml));
}

TEST(BAM_DataSetIO, can_absolutize_resource_paths)
{
    DataSet dataset;
    ReferenceSet referenceDataset;

    dataset.ExternalResources().Add(
        ExternalResource{"PacBio.SubreadFile.SubreadBamFile", "test.fa"});
    referenceDataset.ExternalResources().Add(
        ExternalResource{"PacBio.SubreadFile.SubreadBamFile", "test.fa"});

    const std::string expectedGenericFn{dataset.Path() + "/test.fa"};
    const std::string expectedReferenceFn{referenceDataset.Path() + "/test.fa"};

    std::ostringstream out;
    dataset.SaveToStream(out);
    const std::string genericDatasetXml{out.str()};
    EXPECT_NE(genericDatasetXml.find(expectedGenericFn), std::string::npos);

    out.str("");
    referenceDataset.SaveToStream(out);
    const std::string referenceDatasetXml{out.str()};
    EXPECT_NE(referenceDatasetXml.find(expectedReferenceFn), std::string::npos);
}

TEST(BAM_DataSetIO, ampersands_are_escaped_in_ouput)
{
    SubreadSet ds;
    internal::DataSetElement e{"Description", XsdType::COLLECTION_METADATA};
    e.Text("Transfer location for R&D");
    ds.AddChild(e);

    std::ostringstream out;
    ds.SaveToStream(out);
    EXPECT_TRUE(out.str().find("R&amp;D") != std::string::npos);
}

TEST(BAM_DataSetIO, can_write_relative_paths_in_denovo_datasets)
{
    const std::string file1{"file1.bam"};
    const std::string file2{"subdir/file2.bam"};

    const std::string relativeResource1{"ResourceId=\"file1.bam\""};
    const std::string relativeResource2{"ResourceId=\"subdir/file2.bam\""};

    DataSet dataset;
    dataset.ExternalResources().Add(ExternalResource{"PacBio.SubreadFile.SubreadBamFile", file1});
    dataset.ExternalResources().Add(ExternalResource{"PacBio.SubreadFile.SubreadBamFile", file2});

    std::ostringstream out;
    dataset.SaveToStream(out);

    // contains the expected files, but has appended an absolute path
    const std::string resolvedXml = out.str();
    EXPECT_NE(std::string::npos, resolvedXml.find(file1));
    EXPECT_NE(std::string::npos, resolvedXml.find(file2));
    EXPECT_EQ(std::string::npos, resolvedXml.find(relativeResource1));
    EXPECT_EQ(std::string::npos, resolvedXml.find(relativeResource2));

    out.str("");
    dataset.SaveToStream(out, DataSetPathMode::ALLOW_RELATIVE);

    // contains just the verbatim, relative file path
    const std::string unresolvedXml{out.str()};
    EXPECT_NE(std::string::npos, unresolvedXml.find(file1));
    EXPECT_NE(std::string::npos, unresolvedXml.find(file2));
    EXPECT_NE(std::string::npos, unresolvedXml.find(relativeResource1));
    EXPECT_NE(std::string::npos, unresolvedXml.find(relativeResource2));
}

TEST(BAM_DataSetIO, can_fetch_all_fasta_files)
{
    const DataSet ds{PbbamTestsConfig::Data_Dir + "/fastx/fasta_extension.referenceset.xml"};
    const auto fastaFiles = ds.FastaFiles();
    EXPECT_EQ(3, fastaFiles.size());
    for (const auto& fn : fastaFiles) {
        EXPECT_TRUE(boost::iends_with(fn, "fasta") || boost::iends_with(fn, "fa") ||
                    boost::iends_with(fn, "fsa"));
    }
}

TEST(BAM_DataSetIO, can_read_write_supplemental_resources)
{
    DataSet ds1{PbbamTestsConfig::Data_Dir +
                "/dataset/supplemental_resource1.consensusreadset.xml"};
    const DataSet ds2{PbbamTestsConfig::Data_Dir +
                      "/dataset/supplemental_resource2.consensusreadset.xml"};

    ASSERT_EQ(1, ds1.SupplementalResources().Size());
    EXPECT_EQ("report.txt", ds1.SupplementalResources()[0].ResourceId());
    ASSERT_EQ(1, ds2.SupplementalResources().Size());
    EXPECT_EQ("report2.txt", ds2.SupplementalResources()[0].ResourceId());

    ds1 += ds2;
    ASSERT_EQ(2, ds1.SupplementalResources().Size());
    EXPECT_EQ("report.txt", ds1.SupplementalResources()[0].ResourceId());
    EXPECT_EQ("report2.txt", ds1.SupplementalResources()[1].ResourceId());

    std::ostringstream out;
    ds1.SaveToStream(out);

    const DataSet ds3 = DataSet::FromXml(out.str());
    ASSERT_EQ(2, ds3.SupplementalResources().Size());
    EXPECT_EQ("report.txt", ds3.SupplementalResources()[0].ResourceId());
    EXPECT_EQ("report2.txt", ds3.SupplementalResources()[1].ResourceId());
}

TEST(Bam_DataSetIO, can_merge_from_various_supplemental_resource_counts)
{
    DataSet dataset{PbbamTestsConfig::Data_Dir +
                    "/dataset/supplemental_resource1.consensusreadset.xml"};
    EXPECT_EQ(1, dataset.SupplementalResources().Size());

    const DataSet singleResoureDataset{PbbamTestsConfig::Data_Dir +
                                       "/dataset/supplemental_resource2.consensusreadset.xml"};
    EXPECT_EQ(1, singleResoureDataset.SupplementalResources().Size());

    const DataSet noResourceDataset{PbbamTestsConfig::Data_Dir +
                                    "/dataset/supplemental_resource_empty.consensusreadset.xml"};
    EXPECT_EQ(0, noResourceDataset.SupplementalResources().Size());

    const DataSet mutlipleResourceDataset{
        PbbamTestsConfig::Data_Dir +
        "/dataset/supplemental_resource_multiple.consensusreadset.xml"};
    EXPECT_EQ(3, mutlipleResourceDataset.SupplementalResources().Size());

    dataset += singleResoureDataset;
    EXPECT_EQ(2, dataset.SupplementalResources().Size());

    dataset += noResourceDataset;
    EXPECT_EQ(2, dataset.SupplementalResources().Size());

    dataset += mutlipleResourceDataset;
    EXPECT_EQ(5, dataset.SupplementalResources().Size());
}

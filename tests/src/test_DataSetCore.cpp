// Author: Derek Barnett

#include <cstddef>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/DataSet.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace DataSetCoreTests {

static inline DataSet CreateDataSet()
{
    DataSet d;
    d.Name("foo");
    return d;
}

}  // namespace DataSetCoreTests

TEST(DataSetCoreTest, XmlNameParts)
{
    internal::XmlName name("ns:node_name");
    EXPECT_EQ(boost::string_ref("ns"), name.Prefix());
    EXPECT_EQ(boost::string_ref("node_name"), name.LocalName());
    EXPECT_EQ(boost::string_ref("ns:node_name"), name.QualifiedName());

    internal::XmlName bareName("node_name");
    EXPECT_EQ(boost::string_ref(""), bareName.Prefix());
    EXPECT_EQ(boost::string_ref("node_name"), bareName.LocalName());
    EXPECT_EQ(boost::string_ref("node_name"), bareName.QualifiedName());

    internal::XmlName leadingColon(":node_name");
    EXPECT_EQ(boost::string_ref(""), leadingColon.Prefix());
    EXPECT_EQ(boost::string_ref(":node_name"), leadingColon.LocalName());
    EXPECT_EQ(boost::string_ref(":node_name"), leadingColon.QualifiedName());
}

TEST(DataSetCoreTest, DefaultsOk)
{
    DataSet dataset;
    EXPECT_EQ(DataSet::GENERIC, dataset.Type());
    EXPECT_FALSE(dataset.CreatedAt().empty());
    EXPECT_FALSE(dataset.MetaType().empty());
    EXPECT_FALSE(dataset.TimeStampedName().empty());
    EXPECT_FALSE(dataset.UniqueId().empty());
    EXPECT_FALSE(dataset.Version().empty());

    EXPECT_EQ(0, dataset.TimeStampedName().find("pacbio_dataset_"));

    EXPECT_TRUE(dataset.Format().empty());
    EXPECT_TRUE(dataset.ModifiedAt().empty());
    EXPECT_TRUE(dataset.Name().empty());
    EXPECT_TRUE(dataset.ResourceId().empty());
    EXPECT_TRUE(dataset.Tags().empty());
    EXPECT_EQ(0, dataset.ExternalResources().Size());
    EXPECT_EQ(0, dataset.Filters().Size());
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    EXPECT_EQ(std::string{"3.0.1"}, dataset.Version());
}

TEST(DataSetCoreTest, TimeStampedNamesOk)
{
    DataSet dataset;
    AlignmentSet alignmentSet;
    BarcodeSet barcodeSet;
    ContigSet contigSet;
    ConsensusAlignmentSet consensusAlignmentSet;
    ConsensusReadSet consensusReadSet;
    HdfSubreadSet hdfSubreadSet;
    ReferenceSet referenceSet;
    SubreadSet subreadSet;
    TranscriptSet transcriptSet;

    EXPECT_EQ(0, dataset.TimeStampedName().find("pacbio_dataset_dataset-"));
    EXPECT_EQ(0, alignmentSet.TimeStampedName().find("pacbio_dataset_alignmentset-"));
    EXPECT_EQ(0, barcodeSet.TimeStampedName().find("pacbio_dataset_barcodeset-"));
    EXPECT_EQ(0, contigSet.TimeStampedName().find("pacbio_dataset_contigset-"));
    EXPECT_EQ(
        0, consensusAlignmentSet.TimeStampedName().find("pacbio_dataset_consensusalignmentset-"));
    EXPECT_EQ(0, consensusReadSet.TimeStampedName().find("pacbio_dataset_consensusreadset-"));
    EXPECT_EQ(0, hdfSubreadSet.TimeStampedName().find("pacbio_dataset_hdfsubreadset-"));
    EXPECT_EQ(0, referenceSet.TimeStampedName().find("pacbio_dataset_referenceset-"));
    EXPECT_EQ(0, subreadSet.TimeStampedName().find("pacbio_dataset_subreadset-"));
    EXPECT_EQ(0, transcriptSet.TimeStampedName().find("pacbio_dataset_transcriptset-"));
}

TEST(DataSetCoreTest, BasicGettersSettersOk)
{
    DataSet dataset;
    dataset.CreatedAt("now");
    dataset.Format("format");
    dataset.MetaType("meta");
    dataset.ModifiedAt("later");
    dataset.Name("foo");
    dataset.ResourceId("path/to/file");
    dataset.Tags("tag tag");
    dataset.TimeStampedName("now:30");
    dataset.UniqueId("uuid");
    dataset.Version("0.0.0");

    EXPECT_EQ(std::string("now"), dataset.CreatedAt());
    EXPECT_EQ(std::string("format"), dataset.Format());
    EXPECT_EQ(std::string("meta"), dataset.MetaType());
    EXPECT_EQ(std::string("later"), dataset.ModifiedAt());
    EXPECT_EQ(std::string("foo"), dataset.Name());
    EXPECT_EQ(std::string("path/to/file"), dataset.ResourceId());
    EXPECT_EQ(std::string("tag tag"), dataset.Tags());
    EXPECT_EQ(std::string("now:30"), dataset.TimeStampedName());
    EXPECT_EQ(std::string("uuid"), dataset.UniqueId());
    EXPECT_EQ(std::string("0.0.0"), dataset.Version());
}

TEST(DataSetCoreTest, CopyOk)
{
    DataSet d1;
    d1.Name("foo");

    // copy ctor
    DataSet d2(d1);
    EXPECT_EQ(std::string("foo"), d2.Name());

    // copy assignment
    DataSet d3;
    d3 = d1;
    EXPECT_EQ(std::string("foo"), d3.Name());
}

TEST(DataSetCoreTest, MoveOk)
{
    DataSet d1;
    d1.Name("foo");

// move ctor
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    DataSet d2(std::move(DataSetCoreTests::CreateDataSet()));
#ifdef __clang__
#pragma clang diagnostic pop
#endif
    EXPECT_EQ(std::string("foo"), d2.Name());

    // move assignment
    DataSet d3;
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    d3 = std::move(DataSetCoreTests::CreateDataSet());
#ifdef __clang__
#pragma clang diagnostic pop
#endif
    EXPECT_EQ(std::string("foo"), d3.Name());
}

TEST(DataSetCoreTest, AddExternalResources)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.ExternalResources().Size());

    ExternalResource resource1("metatype", "id");
    resource1.Name("file1");

    ExternalResource resource2("metatype", "id2");
    resource2.Name("file2");

    dataset.ExternalResources().Add(resource1);
    dataset.ExternalResources().Add(resource2);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // disallow duplicates (checking on ResourceId)
    ExternalResource duplicateResource("metatype", "id");
    dataset.ExternalResources().Add(duplicateResource);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // direct access
    const ExternalResources& resources = dataset.ExternalResources();
    EXPECT_EQ(std::string("file1"), resources[0].Name());
    EXPECT_EQ(std::string("file2"), resources[1].Name());

    // iterable
    size_t i = 0;
    for (auto r : resources) {
        if (i == 0)
            EXPECT_EQ(std::string("file1"), r.Name());
        else
            EXPECT_EQ(std::string("file2"), r.Name());
        ++i;
    }
}

TEST(DataSetCoreTest, EditExternalResources)
{
    DataSet dataset;

    ExternalResource resource("metatype", "id");
    resource.Name("file1");
    dataset.ExternalResources().Add(resource);

    resource.Name("file2").ResourceId("id2");
    dataset.ExternalResources().Add(resource);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // edit
    dataset.ExternalResources()[0].Name("some new name");
    EXPECT_EQ(std::string("some new name"), dataset.ExternalResources()[0].Name());
    EXPECT_EQ(std::string("file2"), dataset.ExternalResources()[1].Name());
}

TEST(DataSetCoreTest, NestedExternalResources)
{
    ExternalResource resource("metatype", "filename");
    resource.ExternalResources().Add(ExternalResource("metatype.child", "filename.child"));
    resource.ExternalResources().Add(ExternalResource("metatype.child2", "filename.child2"));

    const ExternalResources& childResources = resource.ExternalResources();
    EXPECT_EQ(2, childResources.Size());
    EXPECT_EQ(std::string("metatype.child"), childResources[0].MetaType());
    EXPECT_EQ(std::string("metatype.child2"), childResources[1].MetaType());
    EXPECT_EQ(std::string("filename.child"), childResources[0].ResourceId());
    EXPECT_EQ(std::string("filename.child2"), childResources[1].ResourceId());
}

TEST(DataSetCoreTest, AddFilters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.Filters().Size());

    Filter filter;
    filter.Properties().Add(Property("rq", "0.85", ">"));
    filter.Properties().Add(Property("RNAME", "chr1", "=="));
    EXPECT_EQ(2, filter.Properties().Size());

    Filter filter2;
    filter2.Properties().Add(Property("rq", "0.50", ">="));
    filter2.Properties().Add(Property("RNAME", "chr2", "!="));
    EXPECT_EQ(2, filter2.Properties().Size());

    dataset.Filters().Add(filter);
    dataset.Filters().Add(filter2);

    const Filters& filters = dataset.Filters();
    EXPECT_EQ(2, filters.Size());
    EXPECT_EQ(2, filters[0].Properties().Size());
    EXPECT_EQ(2, filters[1].Properties().Size());

    // direct access
    const Property& p0 = filters[0].Properties()[0];
    EXPECT_EQ(std::string("rq"), p0.Name());
    EXPECT_EQ(std::string("0.85"), p0.Value());
    EXPECT_EQ(std::string(">"), p0.Operator());

    const Property& p1 = filters[0].Properties()[1];
    EXPECT_EQ(std::string("RNAME"), p1.Name());
    EXPECT_EQ(std::string("chr1"), p1.Value());
    EXPECT_EQ(std::string("=="), p1.Operator());

    const Property& p2 = filters[1].Properties()[0];
    EXPECT_EQ(std::string("rq"), p2.Name());
    EXPECT_EQ(std::string("0.50"), p2.Value());
    EXPECT_EQ(std::string(">="), p2.Operator());

    const Property& p3 = filters[1].Properties()[1];
    EXPECT_EQ(std::string("RNAME"), p3.Name());
    EXPECT_EQ(std::string("chr2"), p3.Value());
    EXPECT_EQ(std::string("!="), p3.Operator());

    // iteratable
    size_t i = 0;
    size_t j = 0;
    for (const Filter& f : filters) {
        if (i == 0) {
            const Properties& properties = f.Properties();
            for (const Property& p : properties) {
                if (j == 0) {
                    EXPECT_EQ(std::string("rq"), p.Name());
                    EXPECT_EQ(std::string("0.85"), p.Value());
                    EXPECT_EQ(std::string(">"), p.Operator());
                } else {
                    EXPECT_EQ(std::string("RNAME"), p.Name());
                    EXPECT_EQ(std::string("chr1"), p.Value());
                    EXPECT_EQ(std::string("=="), p.Operator());
                }
                ++j;
            }
        } else {
            const Properties& properties = f.Properties();
            for (const Property& p : properties) {
                if (j == 0) {
                    EXPECT_EQ(std::string("rq"), p.Name());
                    EXPECT_EQ(std::string("0.50"), p.Value());
                    EXPECT_EQ(std::string(">="), p.Operator());
                } else {
                    EXPECT_EQ(std::string("RNAME"), p.Name());
                    EXPECT_EQ(std::string("chr2"), p.Value());
                    EXPECT_EQ(std::string("!="), p.Operator());
                }
                ++j;
            }
        }
        ++i;
        j = 0;
    }
}

TEST(DataSetCoreTest, EditFilters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.Filters().Size());

    Filter filter;
    filter.Properties().Add(Property("rq", "0.85", ">"));
    filter.Properties().Add(Property("RNAME", "chr1", "=="));
    EXPECT_EQ(2, filter.Properties().Size());

    Filter filter2;
    filter2.Properties().Add(Property("rq", "0.50", ">="));
    filter2.Properties().Add(Property("RNAME", "chr2", "!="));
    EXPECT_EQ(2, filter2.Properties().Size());

    dataset.Filters().Add(filter);
    dataset.Filters().Add(filter2);
    EXPECT_EQ(2, dataset.Filters().Size());
    EXPECT_EQ(2, dataset.Filters()[0].Properties().Size());
    EXPECT_EQ(2, dataset.Filters()[1].Properties().Size());

    // edit property in-place
    Property& p = dataset.Filters()[0].Properties()[0];
    p.Name("someNewName");
    p.Value("someNewValue");
    p.Operator("==");

    const Property& p0 = dataset.Filters()[0].Properties()[0];
    EXPECT_EQ(std::string("someNewName"), p0.Name());
    EXPECT_EQ(std::string("someNewValue"), p0.Value());
    EXPECT_EQ(std::string("=="), p0.Operator());

    const Property& p1 = dataset.Filters()[0].Properties()[1];
    EXPECT_EQ(std::string("RNAME"), p1.Name());
    EXPECT_EQ(std::string("chr1"), p1.Value());
    EXPECT_EQ(std::string("=="), p1.Operator());

    const Property& p2 = dataset.Filters()[1].Properties()[0];
    EXPECT_EQ(std::string("rq"), p2.Name());
    EXPECT_EQ(std::string("0.50"), p2.Value());
    EXPECT_EQ(std::string(">="), p2.Operator());

    const Property& p3 = dataset.Filters()[1].Properties()[1];
    EXPECT_EQ(std::string("RNAME"), p3.Name());
    EXPECT_EQ(std::string("chr2"), p3.Value());
    EXPECT_EQ(std::string("!="), p3.Operator());
}

TEST(DataSetCoreTest, AddSubDataSets)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    DataSetBase sub1;
    sub1.Name("subset_1");

    DataSetBase sub2;
    sub2.Name("subset_2");

    dataset.SubDataSets().Add(sub1);
    dataset.SubDataSets().Add(sub2);
    EXPECT_EQ(2, dataset.SubDataSets().Size());

    // direct access
    const SubDataSets& subdatasets = dataset.SubDataSets();
    EXPECT_EQ(std::string("subset_1"), subdatasets[0].Name());
    EXPECT_EQ(std::string("subset_2"), subdatasets[1].Name());

    // iterable
    size_t i = 0;
    for (const DataSetBase& ds : subdatasets) {
        if (i == 0)
            EXPECT_EQ(std::string("subset_1"), ds.Name());
        else
            EXPECT_EQ(std::string("subset_2"), ds.Name());
        ++i;
    }
}

TEST(DataSetCoreTest, EditSubDataSets)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    DataSetBase sub1;
    sub1.Name("subset_1");

    DataSetBase sub2;
    sub2.Name("subset_2");

    dataset.SubDataSets().Add(sub1);
    dataset.SubDataSets().Add(sub2);
    EXPECT_EQ(2, dataset.SubDataSets().Size());

    // edit
    dataset.SubDataSets()[0].Name("subset_1_edited");

    // direct access
    const SubDataSets& subdatasets = dataset.SubDataSets();
    EXPECT_EQ(std::string("subset_1_edited"), subdatasets[0].Name());
    EXPECT_EQ(std::string("subset_2"), subdatasets[1].Name());

    // iterable
    size_t i = 0;
    for (const DataSetBase& ds : subdatasets) {
        if (i == 0)
            EXPECT_EQ(std::string("subset_1_edited"), ds.Name());
        else
            EXPECT_EQ(std::string("subset_2"), ds.Name());
        ++i;
    }
}

TEST(DataSetCoreTest, RemoveExternalResources)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.ExternalResources().Size());

    ExternalResource resource1("metatype", "id");
    resource1.Name("file1");

    ExternalResource resource2("metatype", "id2");
    resource2.Name("file2");

    dataset.ExternalResources().Add(resource1);
    dataset.ExternalResources().Add(resource2);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // remove
    dataset.ExternalResources().Remove(resource1);
    EXPECT_EQ(1, dataset.ExternalResources().Size());

    // direct access
    const ExternalResources& resources = dataset.ExternalResources();
    EXPECT_EQ(std::string("file2"), resources[0].Name());

    // iterable
    size_t i = 0;
    for (auto r : resources) {
        if (i == 0) {
            EXPECT_EQ(std::string("file2"), r.Name());
        }
        ++i;
    }
}

TEST(DataSetCoreTest, RemoveFilters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.Filters().Size());

    Filter filter;
    filter.Properties().Add(Property("rq", "0.85", ">"));
    filter.Properties().Add(Property("RNAME", "chr1", "=="));
    EXPECT_EQ(2, filter.Properties().Size());

    Filter filter2;
    filter2.Properties().Add(Property("rq", "0.50", ">="));
    filter2.Properties().Add(Property("RNAME", "chr2", "!="));
    EXPECT_EQ(2, filter2.Properties().Size());

    dataset.Filters().Add(filter);
    dataset.Filters().Add(filter2);
    EXPECT_EQ(2, dataset.Filters().Size());

    // remove
    dataset.Filters().Remove(filter);
    EXPECT_EQ(1, dataset.Filters().Size());

    const Filters& filters = dataset.Filters();
    EXPECT_EQ(2, filters[0].Properties().Size());
}

TEST(DataSetCoreTest, RemoveSubDataSets)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.SubDataSets().Size());

    DataSetBase sub1;
    sub1.Name("subset_1");

    DataSetBase sub2;
    sub2.Name("subset_2");

    dataset.SubDataSets().Add(sub1);
    dataset.SubDataSets().Add(sub2);
    EXPECT_EQ(2, dataset.SubDataSets().Size());

    // remove
    dataset.SubDataSets().Remove(sub2);
    EXPECT_EQ(1, dataset.SubDataSets().Size());
}

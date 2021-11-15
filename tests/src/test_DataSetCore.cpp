#include <pbbam/DataSet.h>

#include <cstddef>

#include <string>

#include <gtest/gtest.h>

#include "../src/FileUtils.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace DataSetCoreTests {

const std::string subreadsetBioSample{PbbamTestsConfig::Data_Dir +
                                      "/dataset/biosample.subreadset.xml"};

static DataSet CreateDataSet()
{
    DataSet d;
    d.Name("foo");
    return d;
}

}  // namespace DataSetCoreTests

TEST(BAM_DataSetCore, can_parse_xml_name_parts)
{
    internal::XmlName name{"ns:node_name"};
    EXPECT_EQ(boost::string_ref("ns"), name.Prefix());
    EXPECT_EQ(boost::string_ref("node_name"), name.LocalName());
    EXPECT_EQ(boost::string_ref("ns:node_name"), name.QualifiedName());

    internal::XmlName bareName{"node_name"};
    EXPECT_EQ(boost::string_ref(""), bareName.Prefix());
    EXPECT_EQ(boost::string_ref("node_name"), bareName.LocalName());
    EXPECT_EQ(boost::string_ref("node_name"), bareName.QualifiedName());

    internal::XmlName leadingColon{":node_name"};
    EXPECT_EQ(boost::string_ref(""), leadingColon.Prefix());
    EXPECT_EQ(boost::string_ref(":node_name"), leadingColon.LocalName());
    EXPECT_EQ(boost::string_ref(":node_name"), leadingColon.QualifiedName());
}

TEST(BAM_DataSetCore, created_with_correct_defaults)
{
    const DataSet dataset;
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

TEST(BAM_DataSetCore, default_constructed_generates_time_stamped_name)
{
    const DataSet dataset;
    const AlignmentSet alignmentSet;
    const BarcodeSet barcodeSet;
    const ContigSet contigSet;
    const ConsensusAlignmentSet consensusAlignmentSet;
    const ConsensusReadSet consensusReadSet;
    const HdfSubreadSet hdfSubreadSet;
    const ReferenceSet referenceSet;
    const SubreadSet subreadSet;
    const TranscriptSet transcriptSet;

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

TEST(BAM_DataSetCore, can_be_modified_via_setters)
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

    EXPECT_EQ("now", dataset.CreatedAt());
    EXPECT_EQ("format", dataset.Format());
    EXPECT_EQ("meta", dataset.MetaType());
    EXPECT_EQ("later", dataset.ModifiedAt());
    EXPECT_EQ("foo", dataset.Name());
    EXPECT_EQ("path/to/file", dataset.ResourceId());
    EXPECT_EQ("tag tag", dataset.Tags());
    EXPECT_EQ("now:30", dataset.TimeStampedName());
    EXPECT_EQ("uuid", dataset.UniqueId());
    EXPECT_EQ("0.0.0", dataset.Version());
}

TEST(BAM_DataSetCore, can_be_copied)
{
    DataSet d1;
    d1.Name("foo");

    // copy ctor
    const DataSet d2{d1};
    EXPECT_EQ("foo", d2.Name());

    // copy assignment
    DataSet d3;
    d3 = d1;
    EXPECT_EQ("foo", d3.Name());
}

TEST(BAM_DataSetCore, can_be_moved)
{
    DataSet d1;
    d1.Name("foo");

    // move ctor
    DataSet d2{DataSetCoreTests::CreateDataSet()};
    EXPECT_EQ("foo", d2.Name());

    // move assignment
    DataSet d3;
    d3 = DataSetCoreTests::CreateDataSet();
    EXPECT_EQ("foo", d3.Name());
}

TEST(BAM_DataSetCore, can_add_external_resources)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.ExternalResources().Size());

    ExternalResource resource1{"metatype", "id"};
    resource1.Name("file1");

    ExternalResource resource2{"metatype", "id2"};
    resource2.Name("file2");

    dataset.ExternalResources().Add(resource1);
    dataset.ExternalResources().Add(resource2);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // disallow duplicates (checking on ResourceId)
    const ExternalResource duplicateResource{"metatype", "id"};
    dataset.ExternalResources().Add(duplicateResource);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // direct access
    const ExternalResources& resources = dataset.ExternalResources();
    ASSERT_EQ(2, resources.Size());
    EXPECT_EQ("file1", resources[0].Name());
    EXPECT_EQ("file2", resources[1].Name());

    // iterable
    size_t i = 0;
    for (auto r : resources) {
        if (i == 0) {
            EXPECT_EQ("file1", r.Name());
        } else {
            EXPECT_EQ("file2", r.Name());
        }
        ++i;
    }
}

TEST(BAM_DataSetCore, can_edit_external_resources)
{
    DataSet dataset;

    ExternalResource resource{"metatype", "id"};
    resource.Name("file1");
    dataset.ExternalResources().Add(resource);

    resource.Name("file2").ResourceId("id2");
    dataset.ExternalResources().Add(resource);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // edit
    dataset.ExternalResources()[0].Name("some new name");
    EXPECT_EQ("some new name", dataset.ExternalResources()[0].Name());
    EXPECT_EQ("file2", dataset.ExternalResources()[1].Name());
}

TEST(BAM_DataSetCore, can_create_nested_external_resources)
{
    ExternalResource resource{"metatype", "filename"};
    resource.ExternalResources().Add(ExternalResource{"metatype.child", "filename.child"});
    resource.ExternalResources().Add(ExternalResource{"metatype.child2", "filename.child2"});

    const ExternalResources& childResources = resource.ExternalResources();
    ASSERT_EQ(2, childResources.Size());
    EXPECT_EQ("metatype.child", childResources[0].MetaType());
    EXPECT_EQ("metatype.child2", childResources[1].MetaType());
    EXPECT_EQ("filename.child", childResources[0].ResourceId());
    EXPECT_EQ("filename.child2", childResources[1].ResourceId());
}

TEST(BAM_DataSetCore, can_remove_external_resources)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.ExternalResources().Size());

    ExternalResource resource1{"metatype", "id"};
    resource1.Name("file1");

    ExternalResource resource2{"metatype", "id2"};
    resource2.Name("file2");

    dataset.ExternalResources().Add(resource1);
    dataset.ExternalResources().Add(resource2);
    EXPECT_EQ(2, dataset.ExternalResources().Size());

    // remove
    dataset.ExternalResources().Remove(resource1);
    EXPECT_EQ(1, dataset.ExternalResources().Size());

    // direct access
    const ExternalResources& resources = dataset.ExternalResources();
    EXPECT_EQ("file2", resources[0].Name());

    // iterable
    size_t i = 0;
    for (auto r : resources) {
        if (i == 0) {
            EXPECT_EQ("file2", r.Name());
        }
        ++i;
    }
}

TEST(BAM_DataSetCore, can_add_filters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.Filters().Size());

    Filter filter;
    filter.Properties().Add(Property{"rq", "0.85", ">"});
    filter.Properties().Add(Property{"RNAME", "chr1", "=="});
    EXPECT_EQ(2, filter.Properties().Size());

    Filter filter2;
    filter2.Properties().Add(Property{"rq", "0.50", ">="});
    filter2.Properties().Add(Property{"RNAME", "chr2", "!="});
    EXPECT_EQ(2, filter2.Properties().Size());

    dataset.Filters().Add(filter);
    dataset.Filters().Add(filter2);

    const Filters& filters = dataset.Filters();
    EXPECT_EQ(2, filters.Size());
    EXPECT_EQ(2, filters[0].Properties().Size());
    EXPECT_EQ(2, filters[1].Properties().Size());

    // direct access
    const Property& p0 = filters[0].Properties()[0];
    EXPECT_EQ("rq", p0.Name());
    EXPECT_EQ("0.85", p0.Value());
    EXPECT_EQ(">", p0.Operator());

    const Property& p1 = filters[0].Properties()[1];
    EXPECT_EQ("RNAME", p1.Name());
    EXPECT_EQ("chr1", p1.Value());
    EXPECT_EQ("==", p1.Operator());

    const Property& p2 = filters[1].Properties()[0];
    EXPECT_EQ("rq", p2.Name());
    EXPECT_EQ("0.50", p2.Value());
    EXPECT_EQ(">=", p2.Operator());

    const Property& p3 = filters[1].Properties()[1];
    EXPECT_EQ("RNAME", p3.Name());
    EXPECT_EQ("chr2", p3.Value());
    EXPECT_EQ("!=", p3.Operator());

    // iteratable
    size_t i = 0;
    size_t j = 0;
    for (const Filter& f : filters) {
        if (i == 0) {
            const Properties& properties = f.Properties();
            for (const Property& p : properties) {
                if (j == 0) {
                    EXPECT_EQ("rq", p.Name());
                    EXPECT_EQ("0.85", p.Value());
                    EXPECT_EQ(">", p.Operator());
                } else {
                    EXPECT_EQ("RNAME", p.Name());
                    EXPECT_EQ("chr1", p.Value());
                    EXPECT_EQ("==", p.Operator());
                }
                ++j;
            }
        } else {
            const Properties& properties = f.Properties();
            for (const Property& p : properties) {
                if (j == 0) {
                    EXPECT_EQ("rq", p.Name());
                    EXPECT_EQ("0.50", p.Value());
                    EXPECT_EQ(">=", p.Operator());
                } else {
                    EXPECT_EQ("RNAME", p.Name());
                    EXPECT_EQ("chr2", p.Value());
                    EXPECT_EQ("!=", p.Operator());
                }
                ++j;
            }
        }
        ++i;
        j = 0;
    }
}

TEST(BAM_DataSetCore, can_edit_filters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.Filters().Size());

    Filter filter;
    filter.Properties().Add(Property{"rq", "0.85", ">"});
    filter.Properties().Add(Property{"RNAME", "chr1", "=="});
    EXPECT_EQ(2, filter.Properties().Size());

    Filter filter2;
    filter2.Properties().Add(Property{"rq", "0.50", ">="});
    filter2.Properties().Add(Property{"RNAME", "chr2", "!="});
    EXPECT_EQ(2, filter2.Properties().Size());

    dataset.Filters().Add(filter);
    dataset.Filters().Add(filter2);
    ASSERT_EQ(2, dataset.Filters().Size());
    EXPECT_EQ(2, dataset.Filters()[0].Properties().Size());
    EXPECT_EQ(2, dataset.Filters()[1].Properties().Size());

    // edit property in-place
    Property& p = dataset.Filters()[0].Properties()[0];
    p.Name("someNewName");
    p.Value("someNewValue");
    p.Operator("==");

    const Property& p0 = dataset.Filters()[0].Properties()[0];
    EXPECT_EQ("someNewName", p0.Name());
    EXPECT_EQ("someNewValue", p0.Value());
    EXPECT_EQ("==", p0.Operator());

    const Property& p1 = dataset.Filters()[0].Properties()[1];
    EXPECT_EQ("RNAME", p1.Name());
    EXPECT_EQ("chr1", p1.Value());
    EXPECT_EQ("==", p1.Operator());

    const Property& p2 = dataset.Filters()[1].Properties()[0];
    EXPECT_EQ("rq", p2.Name());
    EXPECT_EQ("0.50", p2.Value());
    EXPECT_EQ(">=", p2.Operator());

    const Property& p3 = dataset.Filters()[1].Properties()[1];
    EXPECT_EQ("RNAME", p3.Name());
    EXPECT_EQ("chr2", p3.Value());
    EXPECT_EQ("!=", p3.Operator());
}

TEST(BAM_DataSetCore, can_remove_filters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.Filters().Size());

    Filter filter;
    filter.Properties().Add(Property{"rq", "0.85", ">"});
    filter.Properties().Add(Property{"RNAME", "chr1", "=="});
    EXPECT_EQ(2, filter.Properties().Size());

    Filter filter2;
    filter2.Properties().Add(Property{"rq", "0.50", ">="});
    filter2.Properties().Add(Property{"RNAME", "chr2", "!="});
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

TEST(BAM_DataSetCore, can_add_subdatasets)
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
    EXPECT_EQ("subset_1", subdatasets[0].Name());
    EXPECT_EQ("subset_2", subdatasets[1].Name());

    // iterable
    size_t i = 0;
    for (const DataSetBase& ds : subdatasets) {
        if (i == 0) {
            EXPECT_EQ("subset_1", ds.Name());
        } else {
            EXPECT_EQ("subset_2", ds.Name());
        }
        ++i;
    }
}

TEST(BAM_DataSetCore, can_edit_subdatasets)
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
    EXPECT_EQ("subset_1_edited", subdatasets[0].Name());
    EXPECT_EQ("subset_2", subdatasets[1].Name());

    // iterable
    size_t i = 0;
    for (const DataSetBase& ds : subdatasets) {
        if (i == 0) {
            EXPECT_EQ("subset_1_edited", ds.Name());
        } else {
            EXPECT_EQ("subset_2", ds.Name());
        }
        ++i;
    }
}

TEST(BAM_DataSetCore, can_remove_subdatasets)
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

TEST(BAM_DataSetCore, generates_created_at_attribute)
{
    const DataSet ds;
    const ReferenceSet ref;
    EXPECT_FALSE(ds.CreatedAt().empty());
    EXPECT_FALSE(ref.CreatedAt().empty());
}

TEST(BAM_DataSetCore, can_add_biosamples)
{
    const std::string barcode_1_1{"lbc1--lbc1"};
    const std::string barcode_1_2{"lbc1--lbc2"};
    const std::string barcode_2_1{"lbc2--lbc1"};
    const std::string barcode_2_2{"lbc2--lbc2"};

    BioSample alice{"Alice"};
    alice.DNABarcodes().Add(barcode_1_1);
    alice.DNABarcodes().Add(barcode_1_2);

    EXPECT_EQ("Alice", alice.Name());
    ASSERT_EQ(2, alice.DNABarcodes().Size());
    EXPECT_EQ(barcode_1_1, alice.DNABarcodes()[0].Name());
    EXPECT_EQ(barcode_1_2, alice.DNABarcodes()[1].Name());
    EXPECT_FALSE(alice.DNABarcodes()[0].UniqueId().empty());
    EXPECT_FALSE(alice.DNABarcodes()[1].UniqueId().empty());

    BioSample bob{"Bob"};
    bob.DNABarcodes().Add(barcode_2_1);
    bob.DNABarcodes().Add(DNABarcode{barcode_2_2, "explicit_uuid"});

    EXPECT_EQ("Bob", bob.Name());
    ASSERT_EQ(2, bob.DNABarcodes().Size());
    EXPECT_EQ(barcode_2_1, bob.DNABarcodes()[0].Name());
    EXPECT_EQ(barcode_2_2, bob.DNABarcodes()[1].Name());
    EXPECT_FALSE(bob.DNABarcodes()[0].UniqueId().empty());
    EXPECT_EQ("explicit_uuid", bob.DNABarcodes()[1].UniqueId());

    DataSet dataset;
    DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(0, metadata.BioSamples().Size());

    metadata.BioSamples().Add(alice);
    metadata.BioSamples().Add(bob);

    ASSERT_EQ(2, metadata.BioSamples().Size());
    EXPECT_EQ("Alice", metadata.BioSamples()[0].Name());
    EXPECT_EQ("Bob", metadata.BioSamples()[1].Name());
}

TEST(BAM_DataSetCore, can_load_biosamples_from_xml)
{
    BAM::DataSet ds{DataSetCoreTests::subreadsetBioSample};
    const auto& metadata = ds.Metadata();
    const auto& biosamples = metadata.BioSamples();

    ASSERT_EQ(1, biosamples.Size());
    EXPECT_EQ("test test", biosamples[0].Name());
}

TEST(BAM_DataSetCore, can_fetch_samples)
{
    const std::set<std::string> expected = {"sample1", "sample2"};
    const BAM::DataSet dataset{PbbamTestsConfig::Data_Dir +
                               "/dataset/samples/dataset_sample_test.subreadset.xml"};
    EXPECT_EQ(dataset.BamFilenames().size(), 3);
    EXPECT_EQ(dataset.Samples(), expected);
}

TEST(BAM_DataSetCore, can_add_supplemental_resources)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.SupplementalResources().Size());

    ExternalResource resource1{"metatype", "id"};
    resource1.Name("file1");

    ExternalResource resource2{"metatype", "id2"};
    resource2.Name("file2");

    dataset.SupplementalResources().Add(resource1);
    dataset.SupplementalResources().Add(resource2);
    EXPECT_EQ(2, dataset.SupplementalResources().Size());

    // disallow duplicates (checking on ResourceId)
    const ExternalResource duplicateResource{"metatype", "id"};
    dataset.SupplementalResources().Add(duplicateResource);
    EXPECT_EQ(2, dataset.SupplementalResources().Size());

    // direct access
    const SupplementalResources& resources = dataset.SupplementalResources();
    ASSERT_EQ(2, resources.Size());
    EXPECT_EQ("file1", resources[0].Name());
    EXPECT_EQ("file2", resources[1].Name());

    // iterable
    size_t i = 0;
    for (auto r : resources) {
        if (i == 0) {
            EXPECT_EQ("file1", r.Name());
        } else {
            EXPECT_EQ("file2", r.Name());
        }
        ++i;
    }
}

TEST(BAM_DataSetCore, can_edit_supplemental_resources)
{
    DataSet dataset;

    ExternalResource resource{"metatype", "id"};
    resource.Name("file1");
    dataset.SupplementalResources().Add(resource);

    resource.Name("file2").ResourceId("id2");
    dataset.SupplementalResources().Add(resource);
    EXPECT_EQ(2, dataset.SupplementalResources().Size());

    // edit
    dataset.SupplementalResources()[0].Name("some new name");
    EXPECT_EQ("some new name", dataset.SupplementalResources()[0].Name());
    EXPECT_EQ("file2", dataset.SupplementalResources()[1].Name());
}

TEST(BAM_DataSetCore, can_remove_supplemental_resources)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.SupplementalResources().Size());

    ExternalResource resource1{"metatype", "id"};
    resource1.Name("file1");

    ExternalResource resource2{"metatype", "id2"};
    resource2.Name("file2");

    dataset.SupplementalResources().Add(resource1);
    dataset.SupplementalResources().Add(resource2);
    EXPECT_EQ(2, dataset.SupplementalResources().Size());

    // remove
    dataset.SupplementalResources().Remove(resource1);
    EXPECT_EQ(1, dataset.SupplementalResources().Size());

    // direct access
    const SupplementalResources& resources = dataset.SupplementalResources();
    EXPECT_EQ("file2", resources[0].Name());

    // iterable
    size_t i = 0;
    for (auto r : resources) {
        if (i == 0) {
            EXPECT_EQ("file2", r.Name());
        }
        ++i;
    }
}

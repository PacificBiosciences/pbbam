#include <pbbam/DataSetXsd.h>

#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/DataSet.h>

#include "PbbamTestData.h"

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_DataSetXsd, populates_default_namespaces)
{
    const NamespaceRegistry registry;

    const NamespaceInfo& baseInfo = registry.Namespace(XsdType::BASE_DATA_MODEL);
    const NamespaceInfo& dsInfo   = registry.Namespace(XsdType::DATASETS);
    const NamespaceInfo& defaultInfo = registry.DefaultNamespace();

    EXPECT_EQ(XsdType::DATASETS, registry.DefaultXsd());

    EXPECT_EQ("pbds",   dsInfo.Name());
    EXPECT_EQ("pbbase", baseInfo.Name());
    EXPECT_EQ("pbds",   defaultInfo.Name());

    EXPECT_EQ("http://pacificbiosciences.com/PacBioBaseDataModel.xsd", baseInfo.Uri());
    EXPECT_EQ("http://pacificbiosciences.com/PacBioDatasets.xsd",      dsInfo.Uri());
    EXPECT_EQ("http://pacificbiosciences.com/PacBioDatasets.xsd",      defaultInfo.Uri());
}

TEST(BAM_DataSetXsd, can_reset_default_namespace)
{
    NamespaceRegistry registry;
    registry.SetDefaultXsd(XsdType::DATASETS);

    const NamespaceInfo& defaultInfo = registry.DefaultNamespace();

    EXPECT_EQ(XsdType::DATASETS, registry.DefaultXsd());
    EXPECT_EQ("pbds", defaultInfo.Name());
    EXPECT_EQ("http://pacificbiosciences.com/PacBioDatasets.xsd", defaultInfo.Uri());
}

TEST(BAM_DataSetXsd, can_edit_namespace_registry)
{
    NamespaceRegistry registry;
    registry.Register(XsdType::DATASETS, NamespaceInfo{"custom", "http://custom/uri.xsd"});

    const NamespaceInfo& dsInfo = registry.Namespace(XsdType::DATASETS);

    EXPECT_EQ("custom",                dsInfo.Name());
    EXPECT_EQ("http://custom/uri.xsd", dsInfo.Uri());
}

TEST(BAM_DataSetXsd, edited_registry_reflected_in_output_xml)
{
    DataSet dataset{DataSet::ALIGNMENT};
    dataset.CreatedAt("2015-01-27T09:00:01");
    dataset.MetaType("PacBio.DataSet.AlignmentSet");
    dataset.Name("DataSet_AlignmentSet");
    dataset.Tags("barcode moreTags mapping mytags");
    dataset.TimeStampedName("my_time_stamped_name");
    dataset.UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c");
    dataset.Attribute("xmlns",              "http://pacificbiosciences.com/PacBioDatasets.xsd")
           .Attribute("xmlns:xsi",          "http://www.w3.org/2001/XMLSchema-instance")
           .Attribute("xsi:schemaLocation", "http://pacificbiosciences.com/PacBioDatasets.xsd");

    ExternalResource ext{"Fake.MetaType", "filename"};
    ext.CreatedAt("2015-01-27T09:00:01");
    ext.TimeStampedName("custom_tsn")
       .UniqueId("my_uuid");
    dataset.ExternalResources().Add(ext);

    dataset.Namespaces().Register(XsdType::BASE_DATA_MODEL, NamespaceInfo{"custom", "http://custom/uri.xsd"});

    std::ostringstream s;
    dataset.SaveToStream(s);
    std::string result = s.str();
    EXPECT_NE(result.find("custom:ExternalResource"), std::string::npos);
}

TEST(BAM_DataSetXsd, namespaces_are_propagated_to_child_elements)
{
    { // default namespaces

        DataSet ds;

        // append child elements that do not have a C++ built-in, nor namespace prefix with addition
        DataSetMetadata& metadata = ds.Metadata();
        metadata.AddChild(internal::DataSetElement{"SummaryStats"});
        metadata.AddChild(internal::DataSetElement{"CopyFiles"});
        metadata.AddChild(internal::DataSetElement{"BioSamples"});
        metadata.AddChild(internal::DataSetElement{"AutomationParameters"});

        std::ostringstream s;
        ds.SaveToStream(s);
        const std::string output = s.str();

        // check that default namespace is propagated properly
        EXPECT_TRUE(output.find("pbds:SummaryStats") != std::string::npos);
        EXPECT_TRUE(output.find("pbmeta:CopyFiles") != std::string::npos);
        EXPECT_TRUE(output.find("pbsample:BioSamples") != std::string::npos);
        EXPECT_TRUE(output.find("pbbase:AutomationParameters") != std::string::npos);
    }

    { // custom namespaces

        DataSet ds;

        // setup custom namespaces
        ds.Namespaces().Register(XsdType::BASE_DATA_MODEL,     NamespaceInfo{"custom_base",   "http://custom/base.xsd"});
        ds.Namespaces().Register(XsdType::COLLECTION_METADATA, NamespaceInfo{"custom_meta",   "http://custom/meta.xsd"});
        ds.Namespaces().Register(XsdType::DATASETS,            NamespaceInfo{"custom_ds",     "http://custom/datasets.xsd"});
        ds.Namespaces().Register(XsdType::SAMPLE_INFO,         NamespaceInfo{"custom_sample", "http://custom/base.xsd"});

        // append child elements that do not have a C++ built-in, nor namespace prefix with addition
        DataSetMetadata& metadata = ds.Metadata();
        metadata.AddChild(internal::DataSetElement{"SummaryStats"});
        metadata.AddChild(internal::DataSetElement{"CopyFiles"});
        metadata.AddChild(internal::DataSetElement{"BioSamples"});
        metadata.AddChild(internal::DataSetElement{"AutomationParameters"});

        std::ostringstream s;
        ds.SaveToStream(s);
        const std::string output = s.str();

        // check that custom namespace is propagated properly
        EXPECT_TRUE(output.find("custom_ds:SummaryStats") != std::string::npos);
        EXPECT_TRUE(output.find("custom_meta:CopyFiles") != std::string::npos);
        EXPECT_TRUE(output.find("custom_sample:BioSamples") != std::string::npos);
        EXPECT_TRUE(output.find("custom_base:AutomationParameters") != std::string::npos);
    }
}

// clang-format on

// Author: Derek Barnett

#include <pbbam/CollectionMetadata.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

#include <gtest/gtest.h>

#include <pbbam/DataSet.h>

#include "PbbamTestData.h"

TEST(BAM_CollectionMetadata, throws_on_empty_xml)
{
    EXPECT_THROW(
        { const auto collectionMetadata = PacBio::BAM::CollectionMetadata::FromRawXml(""); },
        std::runtime_error);
}

TEST(BAM_CollectionMetadata, throws_on_invalid_xml)
{
    EXPECT_THROW(
        { const auto collectionMetadata = PacBio::BAM::CollectionMetadata::FromRawXml("bad xml"); },
        std::runtime_error);
}

TEST(BAM_CollectionMetadata, can_create_from_raw_text_and_attach_to_dataset)
{
    // Create CollectionMetadata from raw XML
    const std::string xmlFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                            "/run_metadata/collection_metadata.xml"};
    std::ifstream xmlIn{xmlFn};
    std::ostringstream xmlText;
    xmlText << xmlIn.rdbuf();

    const auto collectionMetadata = PacBio::BAM::CollectionMetadata::FromRawXml(xmlText.str());
    ASSERT_TRUE(collectionMetadata.HasAttribute("InstrumentName"));
    EXPECT_EQ("Sequel-54076", collectionMetadata.Attribute("InstrumentName"));

    // Load up some existing dataset, check collection metdata
    PacBio::BAM::DataSet subreadSet{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                                    "/run_metadata/barcodes.subreadset.xml"};
    EXPECT_EQ("64008", subreadSet.Metadata().CollectionMetadata().Attribute("InstrumentName"));

    // Attach our new CollectionMetadata and check
    subreadSet.Metadata().CollectionMetadata(collectionMetadata);
    EXPECT_EQ("Sequel-54076",
              subreadSet.Metadata().CollectionMetadata().Attribute("InstrumentName"));
}

TEST(BAM_CollectionMetadata, output_correct_biosample)
{
    // Create CollectionMetadata from raw XML
    const std::string xmlFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                            "/run_metadata/collection_metadata.xml"};
    std::ifstream xmlIn{xmlFn};
    std::ostringstream xmlText;
    xmlText << xmlIn.rdbuf();
    ASSERT_TRUE(xmlText.str().find("<BioSample") != std::string::npos);
    const auto collectionMetadata = PacBio::BAM::CollectionMetadata::FromRawXml(xmlText.str());

    // Load up some existing dataset, attach our collection metdata & check output
    PacBio::BAM::DataSet subreadSet{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                                    "/run_metadata/barcodes.subreadset.xml"};
    subreadSet.Metadata().CollectionMetadata(collectionMetadata);

    std::ostringstream xmlOut;
    subreadSet.SaveToStream(xmlOut);
    EXPECT_EQ(xmlOut.str().find("<pbmeta:BioSample"), std::string::npos);
    EXPECT_NE(xmlOut.str().find("<pbsample:BioSample"), std::string::npos);
}

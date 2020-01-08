// Author: Derek Barnett

#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>

#include <pbbam/RunMetadata.h>

#include "PbbamTestData.h"

// clang-format off
TEST(RunMetadataTest, throws_on_invalid_xml)
{
    auto expectThrow = [](const std::string& xml)
    {
        std::istringstream in{xml};
        EXPECT_THROW(
        {
            const auto c = PacBio::BAM::RunMetadata::Collection(in);
        }, std::runtime_error);
    };

    {
        SCOPED_TRACE("Empty XML");
        expectThrow("");
    }
    {
        SCOPED_TRACE("Incorrect root");
        const std::string xml = R"(<Invalid />)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'ExperimentContainer'");
        const std::string xml = R"(
            <PacBioDataModel>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'Runs'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer />
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'Run'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs />
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'Outputs'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs>
                        <Run />
                    </Runs>
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'SubreadSets'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs>
                        <Run>
                            <Outputs />
                        </Run>
                    </Runs>
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'SubreadSet'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs>
                        <Run>
                            <Outputs>
                                <SubreadSets />
                            </Outputs>
                        </Run>
                    </Runs>
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'DataSetMetadata'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs>
                        <Run>
                            <Outputs>
                                <SubreadSets>
                                    <SubreadSet />
                                </SubreadSets>
                            </Outputs>
                        </Run>
                    </Runs>
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'Collections'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs>
                        <Run>
                            <Outputs>
                                <SubreadSets>
                                    <SubreadSet>
                                        <DataSetMetadata />
                                    </SubreadSet>
                                </SubreadSets>
                            </Outputs>
                        </Run>
                    </Runs>
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
    {
        SCOPED_TRACE("Missing 'CollectionMetadata'");
        const std::string xml = R"(
            <PacBioDataModel>
                <ExperimentContainer>
                    <Runs>
                        <Run>
                            <Outputs>
                                <SubreadSets>
                                    <SubreadSet>
                                        <DataSetMetadata>
                                            <Collections />
                                        </DataSetMetadata>
                                    </SubreadSet>
                                </SubreadSets>
                            </Outputs>
                        </Run>
                    </Runs>
                </ExperimentContainer>
            </PacBioDataModel>)";
        expectThrow(xml);
    }
}

// clang-format on

TEST(RunMetadataTest, can_load_single_collection_from_xml_file)
{
    const std::string xmlFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                            "/run_metadata/id.metadata.xml"};

    // -- CollectionMetadata --
    const auto collection = PacBio::BAM::RunMetadata::Collection(xmlFn);
    EXPECT_EQ("Hydrav2-8A-2-Cell2", collection.SubreadSetName());
    ASSERT_TRUE(collection.HasAutomationParameters());
    ASSERT_TRUE(collection.HasBindingKit());
    ASSERT_TRUE(collection.HasControlKit());
    ASSERT_TRUE(collection.HasSequencingKitPlate());
    ASSERT_TRUE(collection.HasTemplatePrepKit());

    // -- AutomationParameters --
    const auto& automationParameters = collection.AutomationParameters();

    // built-ins
    ASSERT_TRUE(automationParameters.HasSNRCut());
    EXPECT_EQ(3.75, automationParameters.SNRCut());
    ASSERT_TRUE(automationParameters.HasInsertSize());
    EXPECT_EQ(2600, automationParameters.InsertSize());

    // generic access
    ASSERT_TRUE(automationParameters.HasParameter("MovieLength"));
    EXPECT_EQ("360", automationParameters.GetParameter("MovieLength"));

    size_t count = 0;
    for (const auto& p : automationParameters) {
        ++count;
        ASSERT_TRUE(automationParameters.HasParameter(p.first));
        EXPECT_EQ(p.second, automationParameters.GetParameter(p.first));
    }
    EXPECT_EQ(16, count);

    // -- ControlKit --
    const auto& controlKit = collection.ControlKit();
    ASSERT_TRUE(controlKit.HasLeftAdapter());
    EXPECT_EQ("TAGAGAGAGAAAAGGAGGAGGAGGCAACAACAACAACTCTCTCTA", controlKit.LeftAdapter());

    ASSERT_TRUE(controlKit.HasRightAdapter());
    EXPECT_EQ("TAGAGAGAGAAAAGGAGGAGGAGGCAACAACAACAACTCTCTCTA", controlKit.RightAdapter());

    ASSERT_TRUE(controlKit.HasSequence());
    EXPECT_TRUE(controlKit.Sequence().find("TGTCTAGGTCATCTCAACGTAGCTTTGACATATAAC") == 0);

    // -- SequencingKitPlate --
    const auto& sequencingKitPlate = collection.SequencingKitPlate();
    ASSERT_TRUE(sequencingKitPlate.HasPartNumber());
    EXPECT_EQ("101-427-800", sequencingKitPlate.PartNumber());

    // -- TemplatePrepKit --
    const auto& templatePrepKit = collection.TemplatePrepKit();
    ASSERT_TRUE(templatePrepKit.HasPartNumber());
    EXPECT_EQ("100-938-900", templatePrepKit.PartNumber());

    ASSERT_TRUE(templatePrepKit.HasLeftAdaptorSequence());
    EXPECT_EQ("ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
              templatePrepKit.LeftAdaptorSequence());

    ASSERT_TRUE(templatePrepKit.HasRightAdaptorSequence());
    EXPECT_EQ("ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT",
              templatePrepKit.RightAdaptorSequence());

    ASSERT_TRUE(templatePrepKit.HasLeftPrimerSequence());
    EXPECT_EQ("aacggaggaggagga", templatePrepKit.LeftPrimerSequence());

    ASSERT_TRUE(templatePrepKit.HasRightPrimerSequence());
    EXPECT_EQ("aacggaggaggagga", templatePrepKit.RightPrimerSequence());
}

TEST(RunMetadataTest, can_load_multiple_collections_from_xml_file)
{
    const std::string xmlFn{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                            "/run_metadata/id.run.metadata.xml"};

    const auto collections = PacBio::BAM::RunMetadata::Collections(xmlFn);
    ASSERT_EQ(2, collections.size());
    EXPECT_TRUE(collections.find("Hydrav2-8A-1-Cell1") != collections.cend());
    EXPECT_TRUE(collections.find("Hydrav2-8A-2-Cell2") != collections.cend());
}

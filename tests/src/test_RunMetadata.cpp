// Author: Derek Barnett

#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>

#include <pbbam/DataSet.h>
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
    EXPECT_TRUE(collection.HasAutomationParameters());
    EXPECT_TRUE(collection.HasBindingKit());
    EXPECT_TRUE(collection.HasControlKit());
    EXPECT_TRUE(collection.HasSequencingKitPlate());
    EXPECT_TRUE(collection.HasTemplatePrepKit());

    // -- AutomationParameters --
    const auto& automationParameters = collection.AutomationParameters();
    ASSERT_TRUE(automationParameters.HasParameter("MovieLength"));
    ASSERT_TRUE(automationParameters.HasSNRCut());
    EXPECT_EQ(3.75, automationParameters.SNRCut());
    ASSERT_TRUE(automationParameters.HasInsertSize());
    EXPECT_EQ(2600, automationParameters.InsertSize());

    // generic parameter access
    ASSERT_TRUE(automationParameters.HasParameter("MovieLength"));
    EXPECT_EQ("360", automationParameters.GetParameter("MovieLength"));

    // iterable parameters
    size_t count = 0;
    for (const auto& p : automationParameters) {
        ++count;
        ASSERT_TRUE(automationParameters.HasParameter(p.Name()));
        EXPECT_EQ(p.Value(), automationParameters.GetParameter(p.Name()));
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

TEST(RunMetadataTest, can_attach_edited_metadata_to_subreadset)
{
    // load run metadata
    const std::string metadataXml{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                                  "/run_metadata/id.metadata.xml"};

    auto c = PacBio::BAM::RunMetadata::Collection(metadataXml);

    // do some edits
    PacBio::BAM::ControlKit& ck = c.ControlKit();
    ck.LeftAdapter("GATTACA");
    ck.RightAdapter("GATTACA");
    ck.Sequence("AACCGGTT");
    ASSERT_EQ("GATTACA", ck.LeftAdapter());
    ASSERT_EQ("GATTACA", ck.RightAdapter());
    ASSERT_EQ("AACCGGTT", ck.Sequence());

    PacBio::BAM::AutomationParameters& parameters = c.AutomationParameters();
    parameters.InsertSize(10000);

    // load subreadset & attach new run metadata
    const std::string originalSubreadsetXml{PacBio::BAM::PbbamTestsConfig::Data_Dir +
                                            "/run_metadata/id.subreadset.xml"};
    PacBio::BAM::DataSet subreadSet{originalSubreadsetXml};
    auto& metadata = subreadSet.Metadata();
    metadata.CollectionMetadata(c);

    // print new dataset contents
    std::ostringstream out;
    subreadSet.SaveToStream(out);

    // check for edits in the new dataset output
    const std::string& output = out.str();

    const std::string adapterSeq =
        R"(&gt;left_adapter\nGATTACA\n&gt;right_adapter\nGATTACA\n&gt;custom_sequence\nAACCGGTT)";
    EXPECT_TRUE(output.find(adapterSeq) != std::string::npos);

    const std::string insertSize =
        R"(<AutomationParameter Name="InsertSize" SimpleValue="10000" ValueDataType="Int32" />)";
    EXPECT_TRUE(output.find(insertSize) != std::string::npos);
}

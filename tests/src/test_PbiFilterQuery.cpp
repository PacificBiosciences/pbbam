// Author: Derek Barnett

#include <pbbam/PbiFilterQuery.h>

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_PbiFilterQuery, can_perform_normal_filtered_queries)
{
    const auto bamFile = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};

    {
        PbiFilterQuery query(PbiQueryLengthFilter{500, Compare::GREATER_THAN_EQUAL}, bamFile);
        EXPECT_EQ(3, query.NumReads());
        EXPECT_EQ(3, std::distance(query.begin(), query.end()));
    }
    {
        // all records aligned to reverse strand && pos >= 9200
        const std::string queryName{
            "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"};
        const Strand strand = Strand::REVERSE;
        const Position start = 9200;
        const auto filter =
            PbiFilter::Intersection({PbiAlignedStrandFilter{strand},
                                     PbiReferenceStartFilter{start, Compare::GREATER_THAN_EQUAL}});

        PbiFilterQuery query(filter, bamFile);
        EXPECT_EQ(1, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(), [&](const BamRecord& r) {
            return r.AlignedStrand() == strand && r.ReferenceStart() >= start &&
                   r.FullName() == queryName;
        });
        EXPECT_EQ(1, count);
    }
    {
        // all records aligned to forward strand && pos >= 9200
        const std::string queryName{
            "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2114_2531"};
        const Strand strand = Strand::FORWARD;
        const Position start = 9200;
        const auto filter =
            PbiFilter::Intersection({PbiAlignedStrandFilter{Strand::FORWARD},
                                     PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterQuery query(filter, bamFile);
        EXPECT_EQ(1, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(), [&](const BamRecord& r) {
            return r.AlignedStrand() == strand && r.ReferenceStart() >= start &&
                   r.FullName() == queryName;
        });
        EXPECT_EQ(1, count);
    }
    {
        // all records from RG ("b89a4406") with numMatches >= 1200
        const auto filter =
            PbiFilter::Intersection({PbiReadGroupFilter{"b89a4406"},
                                     PbiNumMatchesFilter{1200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterQuery query(filter, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(2, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_EQ(std::string("b89a4406"), r.ReadGroupId());
            EXPECT_GE((r.NumMatches()), 1200);
            if (count == 1)
                EXPECT_EQ(
                    std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/"
                                "14743/2579_4055"),
                    r.FullName());
            else {
                if (count == 2) {
                    EXPECT_EQ(
                        std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_"
                                    "X0/14743/4101_5571"),
                        r.FullName());
                }
            }
        }
        EXPECT_EQ(2, count);
    }
}

// clang-format off
TEST(BAM_PbiFilterQuery, can_iterate_zmw_range_from_dataset_input)
{
    const std::string expectedMovieName{"m64004_190414_193017"};

    const DataSet ds(PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml");
    EXPECT_EQ(3, ds.BamFiles().size());

    {  // movie name

        PbiFilterQuery query{PbiMovieNameFilter{expectedMovieName}, ds};
        EXPECT_EQ(1220, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(),
            [&](const BamRecord& r) {
                return r.MovieName() == expectedMovieName;
            });
        EXPECT_EQ(1220, count);
    }

    {  // sequencing chemistries
        std::set<std::string> chems{ds.SequencingChemistries()};
        std::set<std::string> expected{"S/P3-C1/5.0-8M"};
        EXPECT_TRUE(std::equal(chems.begin(), chems.end(), expected.begin()));
    }

    {  // min ZMW

        const int32_t zmw = 54;
        PbiFilterQuery query{PbiZmwFilter{zmw, Compare::GREATER_THAN}, ds};
        EXPECT_EQ(1220, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(),
            [&](const BamRecord& r) {
                return r.HoleNumber() > zmw;
            });
        EXPECT_EQ(1220, count);
    }

    {  // max ZMW

        const int32_t zmw = 1816;
        PbiFilterQuery query{PbiZmwFilter{1816, Compare::LESS_THAN}, ds};
        EXPECT_EQ(150, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(),
            [&](const BamRecord& r) {
                return r.HoleNumber() < zmw;
            });
        EXPECT_EQ(150, count);
    }

    {  // put all together, from DataSet XML

        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        EXPECT_EQ(150, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(),
            [&](const BamRecord& r) {
                return r.HoleNumber() > 54 &&
                       r.HoleNumber() < 1816;
            });
        EXPECT_EQ(150, count);
    }
    {  // empty filter object - should return all records from the same dataset

        PbiFilterQuery query(PbiFilter{}, ds);
        EXPECT_EQ(1220, query.NumReads());
        EXPECT_EQ(1220, std::distance(query.begin(), query.end()));
    }
    {  // no <Filters> element present at all

        const DataSet dsData(PbbamTestsConfig::GeneratedData_Dir +
                             "/chunking_missingfilters.subreadset.xml");

        const PbiFilter filter = PbiFilter::FromDataSet(dsData);
        PbiFilterQuery query(filter, dsData);
        EXPECT_EQ(1220, query.NumReads());
        EXPECT_EQ(1220, std::distance(query.begin(), query.end()));
    }
    {  // <Filters> element contains no child <Filter> elements

        const DataSet dsData(PbbamTestsConfig::GeneratedData_Dir +
                             "/chunking_emptyfilters.subreadset.xml");

        const PbiFilter filter = PbiFilter::FromDataSet(dsData);
        PbiFilterQuery query(filter, dsData);
        EXPECT_EQ(1220, query.NumReads());
        EXPECT_EQ(1220, std::distance(query.begin(), query.end()));
    }
}
// clang-format on

TEST(BAM_PbiFilterQuery, throws_on_missing_pbi_file)
{
    const PbiFilter filter{PbiZmwFilter{31883}};
    const std::string phi29Bam = PbbamTestsConfig::GeneratedData_Dir + "/missing_pbi.bam";
    const std::string hasPbiBam = PbbamTestsConfig::Data_Dir + "/polymerase/production.scraps.bam";

    {  // single file, missing PBI

        EXPECT_THROW(PbiFilterQuery(filter, phi29Bam), std::runtime_error);
    }

    {  // from dataset, all missing PBI

        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        EXPECT_THROW(PbiFilterQuery(filter, ds), std::runtime_error);
    }

    {  // from dataset, mixed PBI presence

        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.ScrapsBamFile", hasPbiBam));
        EXPECT_THROW(PbiFilterQuery(filter, ds), std::runtime_error);
    }
}

TEST(BAM_PbiFilterQuery, can_filter_using_qname_whitelist_file)
{
    const DataSet ds(PbbamTestsConfig::Data_Dir + "/polymerase/qnameFiltered.subreads.dataset.xml");
    const PbiFilter filter = PbiFilter::FromDataSet(ds);
    PbiFilterQuery query(filter, ds);
    EXPECT_EQ(3, query.NumReads());
    EXPECT_EQ(3, std::distance(query.begin(), query.end()));
}

TEST(BAM_PbiFilterQuery, returns_no_records_from_empty_input)
{
    const BamFile file{PbbamTestsConfig::Data_Dir + "/empty.bam"};
    PbiFilterQuery query{PbiFilter{}, file};
    EXPECT_EQ(0, query.NumReads());
    EXPECT_EQ(0, std::distance(query.begin(), query.end()));
}

TEST(BAM_PbiFilterQuery, can_filter_on_barcoded_data)
{
    const BamFile file{PbbamTestsConfig::Data_Dir + "/phi29.bam"};

    // bc_quality == 1
    {
        PbiFilterQuery query{PbiBarcodeQualityFilter{1}, file};
        EXPECT_EQ(120, query.NumReads());
        EXPECT_EQ(120, std::distance(query.begin(), query.end()));
    }

    // bc_quality != 1
    {
        PbiFilterQuery query{PbiBarcodeQualityFilter{1, Compare::NOT_EQUAL}, file};
        EXPECT_EQ(0, query.NumReads());
        EXPECT_EQ(0, std::distance(query.begin(), query.end()));
    }

    // bc_forward == 0
    {
        PbiFilterQuery query{PbiBarcodeForwardFilter{0}, file};
        EXPECT_EQ(40, query.NumReads());
        EXPECT_EQ(40, std::distance(query.begin(), query.end()));
    }

    // bc_forward == [0,2]
    {
        const auto ids = std::vector<int16_t>{0, 2};
        PbiFilterQuery query{PbiBarcodeForwardFilter{ids}, file};
        EXPECT_EQ(80, query.NumReads());
        EXPECT_EQ(80, std::distance(query.begin(), query.end()));
    }

    // bc_reverse != 0
    {
        PbiFilterQuery query{PbiBarcodeReverseFilter{0, Compare::NOT_EQUAL}, file};
        EXPECT_EQ(80, query.NumReads());
        EXPECT_EQ(80, std::distance(query.begin(), query.end()));
    }
}

TEST(BAM_PbiFilterQuery, can_filter_barcodes_from_xml)
{

    const std::string xml_all = R"_XML_(
<?xml version="1.0" encoding="utf-8"?>
<pbds:SubreadSet
   xmlns="http://pacificbiosciences.com/PacBioDatasets.xsd"
   xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
   xmlns:pbbase="http://pacificbiosciences.com/PacBioBaseDataModel.xsd"
   xmlns:pbsample="http://pacificbiosciences.com/PacBioSampleInfo.xsd"
   xmlns:pbmeta="http://pacificbiosciences.com/PacBioCollectionMetadata.xsd"
   xmlns:pbds="http://pacificbiosciences.com/PacBioDatasets.xsd"
   xsi:schemaLocation="http://pacificbiosciences.com/PacBioDataModel.xsd"
   UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe519c"
   TimeStampedName="subreadset_150304_231155"
   MetaType="PacBio.DataSet.SubreadSet"
   Name="DataSet_SubreadSet"
   Tags=""
   Version="3.0.0"
   CreatedAt="2015-01-27T09:00:01">
<pbbase:ExternalResources>
   <pbbase:ExternalResource
       UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5193"
       TimeStampedName="subread_bam_150304_231155"
       MetaType="PacBio.SubreadFile.SubreadBamFile"
       ResourceId="m64004_190414_193017.1.subreads.bam">
       <pbbase:FileIndices>
           <pbbase:FileIndex
               UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5194"
               TimeStampedName="bam_index_150304_231155"
               MetaType="PacBio.Index.PacBioIndex"
               ResourceId="m64004_190414_193017.1.subreads.bam.pbi"/>
       </pbbase:FileIndices>
   </pbbase:ExternalResource>
</pbbase:ExternalResources>
<pbds:Filters>
    <pbds:Filter>
        <pbbase:Properties>
            <pbbase:Property Name="bq" Operator="=" Value="1"/>
        </pbbase:Properties>
    </pbds:Filter>
</pbds:Filters>
</pbds:SubreadSet>
)_XML_";

    const std::string xml_none = R"_XML_(
<?xml version="1.0" encoding="utf-8"?>
<pbds:SubreadSet
   xmlns="http://pacificbiosciences.com/PacBioDatasets.xsd"
   xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
   xmlns:pbbase="http://pacificbiosciences.com/PacBioBaseDataModel.xsd"
   xmlns:pbsample="http://pacificbiosciences.com/PacBioSampleInfo.xsd"
   xmlns:pbmeta="http://pacificbiosciences.com/PacBioCollectionMetadata.xsd"
   xmlns:pbds="http://pacificbiosciences.com/PacBioDatasets.xsd"
   xsi:schemaLocation="http://pacificbiosciences.com/PacBioDataModel.xsd"
   UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe519c"
   TimeStampedName="subreadset_150304_231155"
   MetaType="PacBio.DataSet.SubreadSet"
   Name="DataSet_SubreadSet"
   Tags=""
   Version="3.0.0"
   CreatedAt="2015-01-27T09:00:01">
<pbbase:ExternalResources>
   <pbbase:ExternalResource
       UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5193"
       TimeStampedName="subread_bam_150304_231155"
       MetaType="PacBio.SubreadFile.SubreadBamFile"
       ResourceId="m64004_190414_193017.1.subreads.bam">
       <pbbase:FileIndices>
           <pbbase:FileIndex
               UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5194"
               TimeStampedName="bam_index_150304_231155"
               MetaType="PacBio.Index.PacBioIndex"
               ResourceId="m64004_190414_193017.1.subreads.bam.pbi"/>
       </pbbase:FileIndices>
   </pbbase:ExternalResource>
</pbbase:ExternalResources>
<pbds:Filters>
    <pbds:Filter>
        <pbbase:Properties>
            <pbbase:Property Name="bq" Operator="!=" Value="1"/>
        </pbbase:Properties>
    </pbds:Filter>
</pbds:Filters>
</pbds:SubreadSet>
)_XML_";

    const BamFile file{PbbamTestsConfig::Data_Dir + "/phi29.bam"};

    {  // filter allows all records
        const DataSet ds = DataSet::FromXml(xml_all);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(120, query.NumReads());
        EXPECT_EQ(120, std::distance(query.begin(), query.end()));
    }
    {  // filter allows no records
        const DataSet ds = DataSet::FromXml(xml_none);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(0, query.NumReads());
        EXPECT_EQ(0, std::distance(query.begin(), query.end()));
    }
}

TEST(BAM_PbiFilterQuery, can_filter_read_groups_from_xml)
{
    const BamFile file{PbbamTestsConfig::Data_Dir + "/phi29.bam"};
    const std::string xmlHeader = R"_XML_(
        <?xml version="1.0" encoding="utf-8"?>
        <pbds:SubreadSet
           xmlns="http://pacificbiosciences.com/PacBioDatasets.xsd"
           xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
           xmlns:pbbase="http://pacificbiosciences.com/PacBioBaseDataModel.xsd"
           xmlns:pbsample="http://pacificbiosciences.com/PacBioSampleInfo.xsd"
           xmlns:pbmeta="http://pacificbiosciences.com/PacBioCollectionMetadata.xsd"
           xmlns:pbds="http://pacificbiosciences.com/PacBioDatasets.xsd"
           xsi:schemaLocation="http://pacificbiosciences.com/PacBioDataModel.xsd"
           UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe519c"
           TimeStampedName="subreadset_150304_231155"
           MetaType="PacBio.DataSet.SubreadSet"
           Name="DataSet_SubreadSet"
           Tags=""
           Version="3.0.0"
           CreatedAt="2015-01-27T09:00:01">
        <pbbase:ExternalResources>
           <pbbase:ExternalResource
               UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5193"
               TimeStampedName="subread_bam_150304_231155"
               MetaType="PacBio.SubreadFile.SubreadBamFile"
               ResourceId="phi29.bam">
               <pbbase:FileIndices>
                   <pbbase:FileIndex
                       UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5194"
                       TimeStampedName="bam_index_150304_231155"
                       MetaType="PacBio.Index.PacBioIndex"
                       ResourceId="phi29.bam.pbi"/>
               </pbbase:FileIndices>
           </pbbase:ExternalResource>
        </pbbase:ExternalResources>
        <pbds:Filters>
            <pbds:Filter>
                <pbbase:Properties>)_XML_";

    const std::string xmlFooter = R"_XML_(
                </pbbase:Properties>
            </pbds:Filter>
        </pbds:Filters>
        </pbds:SubreadSet>
        )_XML_";

    {  // equal
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="qid" Operator="==" Value="-1453990154"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(120, query.NumReads());
        EXPECT_EQ(120, std::distance(query.begin(), query.end()));
    }
    {  // not equal
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="qid" Operator="!=" Value="-1453990154"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(0, query.NumReads());
        EXPECT_EQ(0, std::distance(query.begin(), query.end()));
    }
}

TEST(BAM_PbiFilterQuery, can_filter_zmws_from_xml)
{
    const BamFile file{PbbamTestsConfig::Data_Dir + "/phi29.bam"};
    const std::string xmlHeader = R"_XML_(
        <?xml version="1.0" encoding="utf-8"?>
        <pbds:SubreadSet
           xmlns="http://pacificbiosciences.com/PacBioDatasets.xsd"
           xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
           xmlns:pbbase="http://pacificbiosciences.com/PacBioBaseDataModel.xsd"
           xmlns:pbsample="http://pacificbiosciences.com/PacBioSampleInfo.xsd"
           xmlns:pbmeta="http://pacificbiosciences.com/PacBioCollectionMetadata.xsd"
           xmlns:pbds="http://pacificbiosciences.com/PacBioDatasets.xsd"
           xsi:schemaLocation="http://pacificbiosciences.com/PacBioDataModel.xsd"
           UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe519c"
           TimeStampedName="subreadset_150304_231155"
           MetaType="PacBio.DataSet.SubreadSet"
           Name="DataSet_SubreadSet"
           Tags=""
           Version="3.0.0"
           CreatedAt="2015-01-27T09:00:01">
        <pbbase:ExternalResources>
           <pbbase:ExternalResource
               UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5193"
               TimeStampedName="subread_bam_150304_231155"
               MetaType="PacBio.SubreadFile.SubreadBamFile"
               ResourceId="phi29.bam">
               <pbbase:FileIndices>
                   <pbbase:FileIndex
                       UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5194"
                       TimeStampedName="bam_index_150304_231155"
                       MetaType="PacBio.Index.PacBioIndex"
                       ResourceId="phi29.bam.pbi"/>
               </pbbase:FileIndices>
           </pbbase:ExternalResource>
        </pbbase:ExternalResources>
        <pbds:Filters>
            <pbds:Filter>
                <pbbase:Properties>)_XML_";

    const std::string xmlFooter = R"_XML_(
                </pbbase:Properties>
            </pbds:Filter>
        </pbds:Filters>
        </pbds:SubreadSet>
        )_XML_";

    size_t count_30422 = 0;
    size_t count_648 = 0;
    size_t count_17299 = 0;
    size_t count_whitelist = 0;

    {  // 30422
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="30422"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(13, query.NumReads());
        count_30422 = std::distance(query.begin(), query.end());
        EXPECT_EQ(13, count_30422);
    }
    {  // 648
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="648"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(11, query.NumReads());
        count_648 = std::distance(query.begin(), query.end());
        EXPECT_EQ(11, count_648);
    }
    {  // 17299
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="17299"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(4, query.NumReads());
        count_17299 = std::distance(query.begin(), query.end());
        EXPECT_EQ(4, count_17299);
    }
    {  // now check whitelist
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="[30422,648,17299]"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        EXPECT_EQ(28, query.NumReads());
        count_whitelist = std::distance(query.begin(), query.end());
        EXPECT_EQ(28, count_whitelist);

        EXPECT_EQ(count_30422 + count_648 + count_17299, count_whitelist);
    }
}

TEST(BAM_PbiFilterQuery, can_handle_transcript_records)
{
    const std::string transcriptFn = PbbamTestsConfig::Data_Dir + "/transcript.subreads.bam";

    PbiFilterQuery query{PbiFilter{}, transcriptFn};
    for (const auto& b : query)
        EXPECT_TRUE(b.HasHoleNumber());

    {  // zmw whitelist
        const std::vector<int32_t> whitelist = {1, 3};

        std::vector<int32_t> observed;

        PbiFilter filter{PbiZmwFilter{whitelist}};
        PbiFilterQuery queryData{filter, transcriptFn};
        for (const auto& b : queryData) {
            observed.push_back(b.HoleNumber());
        }

        ASSERT_EQ(2, observed.size());
        EXPECT_EQ(1, observed.at(0));
        EXPECT_EQ(3, observed.at(1));
    }
    {  // zmw bounds
        const PbiFilter filter{
            {PbiZmwFilter{2, Compare::GREATER_THAN_EQUAL}, PbiZmwFilter{4, Compare::LESS_THAN}}};

        std::vector<int32_t> observed;

        PbiFilterQuery queryData{filter, transcriptFn};
        for (const auto& b : queryData) {
            observed.push_back(b.HoleNumber());
        }

        ASSERT_EQ(2, observed.size());
        EXPECT_EQ(2, observed.at(0));
        EXPECT_EQ(3, observed.at(1));
    }
    {  // QNAME

        std::vector<int32_t> observed;

        PbiFilter filter{PbiQueryNameFilter{"transcript/2"}};
        PbiFilterQuery queryData{filter, transcriptFn};
        for (const auto& b : queryData) {
            observed.push_back(b.HoleNumber());
        }

        ASSERT_EQ(1, observed.size());
        EXPECT_EQ(2, observed.at(0));
    }
    {  // QNAME whitelist

        const std::vector<std::string> whitelist = {"transcript/1", "transcript/4"};

        std::vector<int32_t> observed;

        PbiFilter filter{PbiQueryNameFilter{whitelist}};
        PbiFilterQuery queryData{filter, transcriptFn};
        for (const auto& b : queryData) {
            observed.push_back(b.HoleNumber());
        }

        ASSERT_EQ(2, observed.size());
        EXPECT_EQ(1, observed.at(0));
        EXPECT_EQ(4, observed.at(1));
    }

    {  // movie name
        std::vector<int32_t> observed;

        PbiFilter filter{PbiMovieNameFilter{"transcript"}};
        PbiFilterQuery queryData{filter, transcriptFn};
        for (const auto& b : queryData) {
            observed.push_back(b.HoleNumber());
        }

        EXPECT_EQ(4, observed.size());
    }

    {  // movie name from DataSet

        const std::string datasetFn = PbbamTestsConfig::Data_Dir + "/transcriptset.xml";

        std::vector<int32_t> observed;

        PacBio::BAM::DataSet ds(datasetFn);
        PacBio::BAM::PbiFilter filter = PacBio::BAM::PbiFilter::FromDataSet(ds);
        PbiFilterQuery queryData{filter, ds};
        for (const auto& b : queryData) {
            observed.push_back(b.HoleNumber());
        }

        EXPECT_EQ(4, observed.size());
    }
}

TEST(BAM_PbiFilterQuery, can_filter_on_barcoded_read_group_id)
{
    const BamFile bamFile{PbbamTestsConfig::Data_Dir + std::string{"/barcoded_read_groups.bam"}};

    {  //  query read group with no barcodes - should catche all, barcoded or not
        const PbiReadGroupFilter filter{"0d7b28fa"};

        PbiFilterQuery query{filter, bamFile};
        EXPECT_EQ(5, query.NumReads());
        EXPECT_EQ(5, std::distance(query.begin(), query.end()));
    }
    {  // query read group with barcode label

        const ReadGroupInfo rg{"0d7b28fa/0--0"};
        const PbiReadGroupFilter filter{rg};

        PbiFilterQuery query{filter, bamFile};
        EXPECT_EQ(1, query.NumReads());
        EXPECT_EQ(1, std::distance(query.begin(), query.end()));
    }
    {  // query multiple read groups with barcode label

        const ReadGroupInfo rg{"0d7b28fa/0--0"};
        const ReadGroupInfo rg1{"0d7b28fa/1--0"};
        const PbiReadGroupFilter filter{std::vector<ReadGroupInfo>{rg, rg1}};

        PbiFilterQuery query{filter, bamFile};
        EXPECT_EQ(2, query.NumReads());
        EXPECT_EQ(2, std::distance(query.begin(), query.end()));
    }
}

// clang-format off
TEST(BAM_PbiFilterQuery, can_reuse_pbi_index_cache)
{
    const DataSet ds(PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml");
    const auto indexCache = MakePbiIndexCache(ds);

    {
        // min ZMW
        const int32_t zmw = 54;

        PbiFilterQuery query{PbiZmwFilter{zmw, Compare::GREATER_THAN}, ds, indexCache};
        EXPECT_EQ(1220, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(),
            [&](const BamRecord& r) {
                return r.HoleNumber() > zmw;
            });
        EXPECT_EQ(1220, count);
    }

    {  // max ZMW
        const int32_t zmw = 1816;

        PbiFilterQuery query{PbiZmwFilter{zmw, Compare::LESS_THAN}, ds, indexCache};
        EXPECT_EQ(150, query.NumReads());
        const auto count = std::count_if(query.begin(), query.end(),
            [&](const BamRecord& r) {
                return r.HoleNumber() < zmw;
            });
        EXPECT_EQ(150, count);
    }
}
// clang-format on

TEST(BAM_PbiFilterQuery, can_filter_on_qname_whitelist_and_blacklist)
{
    const std::string fn{PbbamTestsConfig::Data_Dir + "/dataset/qname_filter.bam"};

    const std::vector<std::string> recordNames{"singleInsertion/0/0_10", "singleInsertion/0/10_20",
                                               "singleInsertion/1/0_10", "singleInsertion/1/10_20"};
    const std::vector<std::string> whitelist{"singleInsertion/0/0_10", "singleInsertion/1/0_10"};
    const std::vector<std::string> blacklist{"singleInsertion/0/10_20", "singleInsertion/1/10_20"};

    {  // sanity check on input
        PbiFilter filter{};
        PbiFilterQuery query{filter, fn};
        EXPECT_EQ(4, query.NumReads());

        size_t i = 0;
        for (const auto& b : query) {
            EXPECT_EQ(recordNames.at(i), b.FullName());
            ++i;
        }
    }
    {
        // whitelist
        PbiFilter filter{PbiQueryNameFilter{whitelist}};
        PbiFilterQuery query{filter, fn};
        EXPECT_EQ(2, query.NumReads());

        size_t i = 0;
        for (const auto& b : query) {
            if (i == 0) {
                EXPECT_EQ(recordNames.at(0), b.FullName());
            } else if (i == 1) {
                EXPECT_EQ(recordNames.at(2), b.FullName());
            }
            ++i;
        }
    }
    {
        // !whitelist
        PbiFilter filter{PbiQueryNameFilter{whitelist, Compare::NOT_CONTAINS}};
        PbiFilterQuery query{filter, fn};
        EXPECT_EQ(2, query.NumReads());

        size_t i = 0;
        for (const auto& b : query) {
            if (i == 0) {
                EXPECT_EQ(recordNames.at(1), b.FullName());
            } else if (i == 1) {
                EXPECT_EQ(recordNames.at(3), b.FullName());
            }
            ++i;
        }
    }
    {
        // blacklist
        PbiFilter filter{PbiQueryNameFilter{blacklist, Compare::NOT_CONTAINS}};
        PbiFilterQuery query{filter, fn};
        EXPECT_EQ(2, query.NumReads());

        size_t i = 0;
        for (const auto& b : query) {
            if (i == 0) {
                EXPECT_EQ(recordNames.at(0), b.FullName());
            } else if (i == 1) {
                EXPECT_EQ(recordNames.at(2), b.FullName());
            }
            ++i;
        }
    }
}

TEST(BAM_PbiFilter, can_filter_by_movie_name_with_barcoded_read_group_ids)
{
    const auto inputXml = PbbamTestsConfig::Data_Dir + "/barcoded_movie_filter/barcoded.xml";
    const DataSet dataset{inputXml};

    {  // no filter, data has 2 reads
        PbiFilterQuery query{PbiFilter{}, dataset};
        EXPECT_EQ(2, query.NumReads());
        EXPECT_EQ(2, std::distance(query.begin(), query.end()));
    }
    {  // dataset filter has 1 movie name (m54006_200116_134114)
        PbiFilterQuery query{PbiFilter::FromDataSet(dataset), dataset};
        EXPECT_EQ(1, query.NumReads());
        EXPECT_EQ(1, std::distance(query.begin(), query.end()));
    }
    {  // use other movie name 'manually'
        PbiFilterQuery query{PbiMovieNameFilter{"m54006_200116_200000"}, dataset};
        EXPECT_EQ(1, query.NumReads());
        EXPECT_EQ(1, std::distance(query.begin(), query.end()));
    }
}

TEST(BAM_PbiFilter, can_filter_subread_records_by_qname)
{
    const std::string fn{
        PbbamTestsConfig::Data_Dir +
        "/chunking/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam"};

    const std::vector<std::string> qnames{"m64004_190414_193017/2865/7276_7872",
                                          "m64004_190414_193017/2865/15855_16411"};

    int count = 0;
    PbiFilterQuery query{PbiQueryNameFilter{qnames}, fn};
    for (const auto& b : query) {
        std::ignore = b;
        ++count;
    }
    EXPECT_EQ(2, count);
}

TEST(BAM_PbiFilter, can_filter_ccs_records_by_qname)
{
    const std::string fn{PbbamTestsConfig::Data_Dir +
                         "/ccs-kinetics-bystrandify/ccs-kinetics-bystrandify-mock-input.2.bam"};

    int count = 0;
    PbiFilterQuery query{PbiQueryNameFilter{"m64011_190228_190319/3/ccs"}, fn};
    for (const auto& b : query) {
        std::ignore = b;
        ++count;
    }
    EXPECT_EQ(1, count);
}

TEST(BAM_PbiFilter, can_filter_transcript_records_by_qname)
{
    const std::string fn{PbbamTestsConfig::Data_Dir + "/transcript.subreads.bam"};

    int count = 0;
    PbiFilterQuery query{PbiQueryNameFilter{"transcript/2"}, fn};
    for (const auto& b : query) {
        std::ignore = b;
        ++count;
    }
    EXPECT_EQ(1, count);
}

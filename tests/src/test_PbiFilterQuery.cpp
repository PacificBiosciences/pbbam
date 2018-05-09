// Author: Derek Barnett

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/PbiFilterQuery.h>
#include <pbbam/Unused.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(PbiFilterQueryTest, QueryOk)
{
    const auto bamFile = BamFile{PbbamTestsConfig::Data_Dir + std::string{"/group/test2.bam"}};

    {
        PbiFilterQuery query(PbiQueryLengthFilter{500, Compare::GREATER_THAN_EQUAL}, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(3, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 500);
        }
        EXPECT_EQ(3, count);
    }
    {
        // all records aligned to reverse strand && pos >= 9200
        const auto filter =
            PbiFilter::Intersection({PbiAlignedStrandFilter{Strand::REVERSE},
                                     PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterQuery query(filter, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(1, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_EQ(Strand::REVERSE, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(
                std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/"
                            "5615_6237"),
                r.FullName());
        }
        EXPECT_EQ(1, count);
    }
    {
        // all records aligned to forward strand && pos >= 9200
        const auto filter =
            PbiFilter::Intersection({PbiAlignedStrandFilter{Strand::FORWARD},
                                     PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}});

        PbiFilterQuery query(filter, bamFile);
        const auto numReads = query.NumReads();
        EXPECT_EQ(1, numReads);

        int count = 0;
        for (const auto& r : query) {
            ++count;
            EXPECT_EQ(Strand::FORWARD, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(
                std::string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/"
                            "2114_2531"),
                r.FullName());
        }
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

TEST(PbiFilterQueryTest, ZmwRangeFromDatasetOk)
{
    const auto expectedMovieName =
        std::string{"m150404_101626_42267_c100807920800000001823174110291514_s1_p0"};

    const DataSet ds(PbbamTestsConfig::Data_Dir + "/chunking/chunking.subreadset.xml");
    EXPECT_EQ(3, ds.BamFiles().size());

    {  // movie name

        PbiFilterQuery query{PbiMovieNameFilter{expectedMovieName}, ds};
        const auto numReads = query.NumReads();
        EXPECT_EQ(1220, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            EXPECT_EQ(expectedMovieName, r.MovieName());
            ++count;
        }
        EXPECT_EQ(1220, count);
    }

    {  // sequencing chemistries
        std::set<std::string> chems{ds.SequencingChemistries()};
        std::set<std::string> expected{"P6-C4"};
        EXPECT_TRUE(equal(chems.begin(), chems.end(), expected.begin()));
    }

    {  // min ZMW

        PbiFilterQuery query{PbiZmwFilter{54, Compare::GREATER_THAN}, ds};
        const auto numReads = query.NumReads();
        EXPECT_EQ(1220, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            EXPECT_GT(r.HoleNumber(), 54);
            ++count;
        }
        EXPECT_EQ(1220, count);
    }

    {  // max ZMW

        PbiFilterQuery query{PbiZmwFilter{1816, Compare::LESS_THAN}, ds};
        const auto numReads = query.NumReads();
        EXPECT_EQ(150, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            EXPECT_LT(r.HoleNumber(), 1816);
            ++count;
        }
        EXPECT_EQ(150, count);
    }

    {  // put all together, from DataSet XML

        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        const auto numReads = query.NumReads();
        EXPECT_EQ(150, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            EXPECT_EQ(expectedMovieName, r.MovieName());
            const auto zmw = r.HoleNumber();
            EXPECT_GT(zmw, 54);
            EXPECT_LT(zmw, 1816);
            ++count;
        }
        EXPECT_EQ(150, count);
    }
    {  // empty filter object - should return all records from the same dataset

        PbiFilterQuery query(PbiFilter{}, ds);
        const auto numReads = query.NumReads();
        EXPECT_EQ(1220, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    {  // no <Filters> element present at all

        const DataSet dsData(PbbamTestsConfig::GeneratedData_Dir +
                             "/chunking_missingfilters.subreadset.xml");
        const PbiFilter filter = PbiFilter::FromDataSet(dsData);
        PbiFilterQuery query(filter, dsData);
        const auto numReads = query.NumReads();
        EXPECT_EQ(1220, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    {  // <Filters> element contains no child <Filter> elements

        const DataSet dsData(PbbamTestsConfig::GeneratedData_Dir +
                             "/chunking_emptyfilters.subreadset.xml");
        const PbiFilter filter = PbiFilter::FromDataSet(dsData);
        PbiFilterQuery query(filter, dsData);
        const auto numReads = query.NumReads();
        EXPECT_EQ(1220, numReads);

        int count = 0;
        for (const BamRecord& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
}

TEST(PbiFilterQueryTest, MissingPbiShouldThrow)
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

TEST(PbiFilterQueryTest, QNameWhitelistFile)
{
    const DataSet ds(PbbamTestsConfig::Data_Dir + "/polymerase/qnameFiltered.subreads.dataset.xml");
    const PbiFilter filter = PbiFilter::FromDataSet(ds);
    PbiFilterQuery query(filter, ds);
    const auto numReads = query.NumReads();
    EXPECT_EQ(3, numReads);

    int count = 0;
    for (const BamRecord& r : query) {
        UNUSED(r);
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(PbiFilterQueryTest, EmptyFiles)
{
    const BamFile file{PbbamTestsConfig::Data_Dir + "/empty.bam"};
    PbiFilterQuery query{PbiFilter{}, file};
    const auto numReads = query.NumReads();
    EXPECT_EQ(0, numReads);

    size_t count = 0;
    for (const auto& r : query) {
        UNUSED(r);
        ++count;
    }
    EXPECT_EQ(0, count);
}

TEST(PbiFilterQueryTest, BarcodeData)
{
    const BamFile file{PbbamTestsConfig::Data_Dir + "/phi29.bam"};

    // bc_quality == 1
    {
        PbiFilterQuery query{PbiBarcodeQualityFilter{1}, file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(120, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(120, count);
    }

    // bc_quality != 1
    {
        PbiFilterQuery query{PbiBarcodeQualityFilter{1, Compare::NOT_EQUAL}, file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(0, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(0, count);
    }

    // bc_forward == 0
    {
        PbiFilterQuery query{PbiBarcodeForwardFilter{0}, file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(40, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(40, count);
    }

    // bc_forward == [0,2]
    {
        const auto ids = std::vector<int16_t>{0, 2};
        PbiFilterQuery query{PbiBarcodeForwardFilter{ids}, file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(80, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(80, count);
    }

    // bc_reverse != 0
    {
        PbiFilterQuery query{PbiBarcodeReverseFilter{0, Compare::NOT_EQUAL}, file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(80, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(80, count);
    }
}

TEST(PbiFilterQueryTest, BarcodeQualityFromXml)
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
       ResourceId="m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam">
       <pbbase:FileIndices>
           <pbbase:FileIndex 
               UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5194" 
               TimeStampedName="bam_index_150304_231155" 
               MetaType="PacBio.Index.PacBioIndex" 
               ResourceId="m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam.pbi"/>
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
       ResourceId="m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam">
       <pbbase:FileIndices>
           <pbbase:FileIndex
               UniqueId="b095d0a3-94b8-4918-b3af-a3f81bbe5194"
               TimeStampedName="bam_index_150304_231155"
               MetaType="PacBio.Index.PacBioIndex"
               ResourceId="m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam.pbi"/>
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
        const auto numReads = query.NumReads();
        EXPECT_EQ(120, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(120, count);
    }
    {  // filter allows no records
        const DataSet ds = DataSet::FromXml(xml_none);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(0, numReads);

        size_t count = 0;
        for (const auto& r : query) {
            UNUSED(r);
            ++count;
        }
        EXPECT_EQ(0, count);
    }
}

TEST(PbiFilterQueryTest, ZmwWhitelistFromXml)
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
        const auto numReads = query.NumReads();
        EXPECT_EQ(13, numReads);

        for (const auto& r : query) {
            UNUSED(r);
            ++count_30422;
        }
        EXPECT_EQ(13, count_30422);
    }
    {  // 648
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="648"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(11, numReads);

        for (const auto& r : query) {
            UNUSED(r);
            ++count_648;
        }
        EXPECT_EQ(11, count_648);
    }
    {  // 17299
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="17299"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(4, numReads);

        for (const auto& r : query) {
            UNUSED(r);
            ++count_17299;
        }
        EXPECT_EQ(4, count_17299);
    }
    {  // now check whitelist
        const std::string xmlProperty =
            R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="[30422,648,17299]"/>\n)_XML_";
        const std::string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query{PbiFilter::FromDataSet(ds), file};
        const auto numReads = query.NumReads();
        EXPECT_EQ(28, numReads);

        for (const auto& r : query) {
            UNUSED(r);
            ++count_whitelist;
        }
        EXPECT_EQ(count_30422 + count_648 + count_17299, count_whitelist);
    }
}

TEST(PbiFilterQueryTest, TranscriptRecords)
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

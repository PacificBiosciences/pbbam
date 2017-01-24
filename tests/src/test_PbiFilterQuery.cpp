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
#include <gtest/gtest.h>
#include <pbbam/PbiFilterQuery.h>
#include <algorithm>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

TEST(PbiFilterQueryTest, QueryOk)
{
    const auto bamFile = BamFile{ tests::Data_Dir + string{ "/group/test2.bam" } };

    {
        int count = 0;
        PbiFilterQuery query( PbiQueryLengthFilter{ 500, Compare::GREATER_THAN_EQUAL}, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_GE((r.QueryEnd() - r.QueryStart()), 500);
        }
        EXPECT_EQ(3, count);
    }
    {
        // all records aligned to reverse strand && pos >= 9200
        const auto filter = PbiFilter::Intersection(
        {
            PbiAlignedStrandFilter{Strand::REVERSE},
            PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}
        });

        int count = 0;
        PbiFilterQuery query(filter, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_EQ(Strand::REVERSE, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/5615_6237"), r.FullName());
        }
        EXPECT_EQ(1, count);
    }
    {
        // all records aligned to forward strand && pos >= 9200
        const auto filter = PbiFilter::Intersection(
        {
            PbiAlignedStrandFilter{Strand::FORWARD},
            PbiReferenceStartFilter{9200, Compare::GREATER_THAN_EQUAL}
        });

        int count = 0;
        PbiFilterQuery query(filter, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_EQ(Strand::FORWARD, r.AlignedStrand());
            EXPECT_GE((r.ReferenceStart()), 9200);
            EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2114_2531"), r.FullName());
        }
        EXPECT_EQ(1, count);
    }
    {
        // all records from RG ("b89a4406") with numMatches >= 1200
        const auto filter = PbiFilter::Intersection(
        {
            PbiReadGroupFilter{"b89a4406"},
            PbiNumMatchesFilter{1200, Compare::GREATER_THAN_EQUAL}
        });

        int count = 0;
        PbiFilterQuery query(filter, bamFile);
        for (const auto& r: query) {
            ++count;
            EXPECT_EQ(string("b89a4406"), r.ReadGroupId());
            EXPECT_GE((r.NumMatches()), 1200);
            if (count == 1)
                EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/2579_4055"), r.FullName());
            else if (count == 2)
                EXPECT_EQ(string("m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/14743/4101_5571"), r.FullName());
        }
        EXPECT_EQ(2, count);
    }
}

TEST(PbiFilterQueryTest, ZmwRangeFromDatasetOk)
{
    const auto expectedMovieName = string{ "m150404_101626_42267_c100807920800000001823174110291514_s1_p0" };

    const DataSet ds(tests::Data_Dir + "/chunking/chunking.subreadset.xml");
    EXPECT_EQ(3, ds.BamFiles().size());

    { // movie name

        int count = 0;
        PbiFilterQuery query{ PbiMovieNameFilter{expectedMovieName}, ds };
        for (const BamRecord& r : query) {
            EXPECT_EQ(expectedMovieName, r.MovieName());
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    
    { // sequencing chemistries
        set<string> chems{ ds.SequencingChemistries() };
        set<string> expected{ "P6-C4" };
        EXPECT_TRUE(equal(chems.begin(), chems.end(), expected.begin()));
    }

    { // min ZMW

        int count = 0;
        PbiFilterQuery query{ PbiZmwFilter{54, Compare::GREATER_THAN}, ds };
        for (const BamRecord& r : query) {
            EXPECT_GT(r.HoleNumber(), 54);
            ++count;
        }
        EXPECT_EQ(1220, count);
    }

    { // max ZMW

        int count = 0;
        PbiFilterQuery query{ PbiZmwFilter{1816, Compare::LESS_THAN}, ds };
        for (const BamRecord& r : query) {
            EXPECT_LT(r.HoleNumber(),1816);
            ++count;
        }
        EXPECT_EQ(150, count);
    }

    { // put all together, from DataSet XML

        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
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
    { // empty filter object - should return all records from the same dataset

        PbiFilterQuery query(PbiFilter{ }, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    { // no <Filters> element present at all

        const DataSet ds(tests::GeneratedData_Dir + "/chunking_missingfilters.subreadset.xml");
        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
    { // <Filters> element contains no child <Filter> elements

        const DataSet ds(tests::GeneratedData_Dir + "/chunking_emptyfilters.subreadset.xml");
        const PbiFilter filter = PbiFilter::FromDataSet(ds);
        PbiFilterQuery query(filter, ds);
        int count = 0;
        for (const BamRecord& r : query) {
            (void)r;
            ++count;
        }
        EXPECT_EQ(1220, count);
    }
}

TEST(PbiFilterQueryTest, MissingPbiShouldThrow)
{
    const PbiFilter filter{ PbiZmwFilter{31883} };
    const string phi29Bam = tests::GeneratedData_Dir + "/missing_pbi.bam";
    const string hasPbiBam = tests::Data_Dir + "/polymerase/production.scraps.bam";

    { // single file, missing PBI

        EXPECT_THROW(PbiFilterQuery(filter, phi29Bam), std::runtime_error);
    }

    { // from dataset, all missing PBI

        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        EXPECT_THROW(PbiFilterQuery(filter, ds), std::runtime_error);
    }

    { // from dataset, mixed PBI presence

        DataSet ds;
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.SubreadBamFile", phi29Bam));
        ds.ExternalResources().Add(ExternalResource("PacBio.SubreadFile.ScrapsBamFile", hasPbiBam));
        EXPECT_THROW(PbiFilterQuery(filter, ds), std::runtime_error);
    }
}

TEST(PbiFilterQueryTest, QNameWhitelistFile)
{
    const DataSet ds(tests::Data_Dir + "/polymerase/qnameFiltered.subreads.dataset.xml");
    const PbiFilter filter = PbiFilter::FromDataSet(ds);
    PbiFilterQuery query(filter, ds);
    int count = 0;
    for (const BamRecord& r : query) {
        (void)r;
        ++count;
    }
    EXPECT_EQ(3, count);
}

TEST(PbiFilterQueryTest, EmptyFiles)
{
    const BamFile file{ tests::Data_Dir + "/empty.bam" };
    PbiFilterQuery query{ PbiFilter{}, file };
    size_t count = 0;
    for (const auto& r : query) {
        (void)r;
        ++count;
    }
    EXPECT_EQ(0, count);
}

TEST(PbiFilterQueryTest, BarcodeData)
{
   const BamFile file{ tests::Data_Dir + "/phi29.bam" };

   // bc_quality == 1
   {
       size_t count = 0;
       PbiFilterQuery query{ PbiBarcodeQualityFilter{1}, file };
       for (const auto& r : query) {
           (void)r;
           ++count;
       }
       EXPECT_EQ(120, count);
   }

   // bc_quality != 1
   {
       size_t count = 0;
       PbiFilterQuery query{ PbiBarcodeQualityFilter{1, Compare::NOT_EQUAL}, file };
       for (const auto& r : query) {
           (void)r;
           ++count;
       }
       EXPECT_EQ(0, count);
   }

   // bc_forward == 0
   {
       size_t count = 0;
       PbiFilterQuery query{ PbiBarcodeForwardFilter{0}, file };
       for (const auto& r : query) {
           (void)r;
           ++count;
       }
       EXPECT_EQ(40, count);
   }

   // bc_forward == [0,2]
   {
       size_t count = 0;
       const auto ids = vector<int16_t>{ 0, 2 };
       PbiFilterQuery query{ PbiBarcodeForwardFilter{ ids }, file };
       for (const auto& r : query) {
           (void)r;
           ++count;
       }
       EXPECT_EQ(80, count);
   }

   // bc_reverse != 0
   {
       size_t count = 0;
       PbiFilterQuery query{ PbiBarcodeReverseFilter{0, Compare::NOT_EQUAL}, file };
       for (const auto& r : query) {
           (void)r;
           ++count;
       }
       EXPECT_EQ(80, count);
   }
}

TEST(PbiFilterQueryTest, BarcodeQualityFromXml)
{

const string xml_all = R"_XML_(
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

const string xml_none = R"_XML_(
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

     const BamFile file{ tests::Data_Dir + "/phi29.bam" }; 

     {   // filter allows all records
         const DataSet ds = DataSet::FromXml(xml_all);
         const PbiFilterQuery query { PbiFilter::FromDataSet(ds), file };
         size_t count = 0;
         for (const auto& r : query) {
             (void)r;
             ++count;
         }
         EXPECT_EQ(120, count);
     }
     {    // filter allows no records
         const DataSet ds = DataSet::FromXml(xml_none);
         const PbiFilterQuery query { PbiFilter::FromDataSet(ds), file };
         size_t count = 0;
         for (const auto& r : query) {
             (void)r;
             ++count;
         }
         EXPECT_EQ(0, count);
    }
}

TEST(PbiFilterQueryTest, ZmwWhitelistFromXml)
{
    const BamFile file{ tests::Data_Dir + "/phi29.bam" };
    const string xmlHeader = R"_XML_(
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

        const string xmlFooter = R"_XML_(
                </pbbase:Properties>
            </pbds:Filter>
        </pbds:Filters>
        </pbds:SubreadSet>
        )_XML_";

    size_t count_30422 = 0;
    size_t count_648 = 0;
    size_t count_17299 = 0;
    size_t count_whitelist = 0;

    {   // 30422
        const string xmlProperty = R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="30422"/>\n)_XML_";
        const string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query { PbiFilter::FromDataSet(ds), file };
        for (const auto& r : query) {
            (void)r;
            ++count_30422;
        }
        EXPECT_EQ(13, count_30422);
    }
    {   // 648
        const string xmlProperty = R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="648"/>\n)_XML_";
        const string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query { PbiFilter::FromDataSet(ds), file };
        for (const auto& r : query) {
            (void)r;
            ++count_648;
        }
        EXPECT_EQ(11, count_648);
    }
    {   // 17299
        const string xmlProperty = R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="17299"/>\n)_XML_";
        const string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query { PbiFilter::FromDataSet(ds), file };
        for (const auto& r : query) {
            (void)r;
            ++count_17299;
        }
        EXPECT_EQ(4, count_17299);
    }
    {   // now check whitelist
        const string xmlProperty = R"_XML_(<pbbase:Property Name="zm" Operator="=" Value="[30422,648,17299]"/>\n)_XML_";
        const string xml = xmlHeader + xmlProperty + xmlFooter;
        const DataSet ds = DataSet::FromXml(xml);
        const PbiFilterQuery query { PbiFilter::FromDataSet(ds), file };
        for (const auto& r : query) {
            (void)r;
            ++count_whitelist;
        }
        EXPECT_EQ(count_30422 + count_648 + count_17299, count_whitelist);
    }
}

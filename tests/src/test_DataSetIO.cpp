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
#include <pbbam/dataset/DataSet.h>
#include <pbbam/internal/DataSetElement.h>
#include <stdexcept>
#include <sstream>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

const string ex2BamFn      = tests::Data_Dir + "/ex2.bam";
const string bamGroupFofn  = tests::Data_Dir + "/test_group_query/group.fofn";

const string ali2XmlFn     = tests::Data_Dir + "/dataset/ali2.xml";
const string subread1XmlFn = tests::Data_Dir + "/dataset/subread_dataset1.xml";
const string exampleXmlFn  = tests::Data_Dir + "/dataset/ex1.xml";

static void TestExampleXml(void);
static void TestSubread1Xml(void);
static void TestAli2Xml(void);

namespace tests {

} // namespace tests

TEST(DataSetIOTest, FromBamFilename)
{
    DataSet dataset(ex2BamFn);

    EXPECT_EQ(1, dataset.ExternalDataReferenceList().Size());
    ExternalDataReference bamRef = dataset.ExternalDataReferenceList()[0];

    EXPECT_EQ(ex2BamFn, bamRef.ResourceId());
}

TEST(DataSetIOTest, FromBamFileObject)
{
    BamFile bamFile(ex2BamFn);
    DataSet dataset(bamFile);

    EXPECT_EQ(1, dataset.ExternalDataReferenceList().Size());
    ExternalDataReference bamRef = dataset.ExternalDataReferenceList()[0];

    EXPECT_EQ(ex2BamFn, bamRef.ResourceId());
}

TEST(DataSetIOTest, FromFofn)
{
    DataSet dataset(bamGroupFofn);

    EXPECT_EQ(3, dataset.ExternalDataReferenceList().Size());
}

TEST(DataSetIOTest, FromXml)
{
    EXPECT_NO_THROW(TestExampleXml());
    EXPECT_NO_THROW(TestSubread1Xml());
    EXPECT_NO_THROW(TestAli2Xml());
}

static void TestAli2Xml(void)
{
    const DataSet dataset(ali2XmlFn);

    // header attributes
    EXPECT_EQ(DataSetType::ALIGNMENTSET,                      dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),                  dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.AlignmentSet"),          dataset.MetaType());
    EXPECT_EQ(string("DataSet_AlignmentSet"),                 dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"),      dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                dataset.Version());

    // external refs
    const ExternalDataReferences& refs = dataset.ExternalDataReferenceList();
    EXPECT_EQ(2, refs.Size());

    ExternalDataReference ref;
    ref = refs[0];
    EXPECT_EQ(string("Points to an example Alignments BAM file."), ref.Description());
    EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"),     ref.MetaType());
    EXPECT_EQ(string("Third Alignments BAM"),               ref.Name());
    EXPECT_EQ(string("file:/mnt/path/to/alignments2.bam"),  ref.ResourceId());
    EXPECT_EQ(string("Example"),                            ref.Tags());

    ref = refs[1];
    EXPECT_EQ(string("Points to another example Alignments BAM file, by relative path."), ref.Description());
    EXPECT_EQ(string("AlignmentFile.AlignmentBamFile"), ref.MetaType());
    EXPECT_EQ(string("Fourth Alignments BAM"),          ref.Name());
    EXPECT_EQ(string("file:./alignments3.bam"),         ref.ResourceId());
    EXPECT_EQ(string("Example"),                        ref.Tags());

    // filters
    EXPECT_EQ(0, dataset.NumFilters());

    // subdatasets
    const SubDataSets& subdatasets = dataset.SubDataSetList();
    EXPECT_EQ(2, subdatasets.Size());

    const SubDataSet& subdataset1 = subdatasets[0];
    EXPECT_EQ(string("HighQuality Read Alignments"),          subdataset1.Name());
    EXPECT_EQ(string("ab95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset1.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                subdataset1.Version());

    const Filters& subFilters1 = subdataset1.FilterList();
    EXPECT_EQ(1, subFilters1.Size());

    const Filter& subFilter1 = subFilters1[0];
    const FilterParameters& subFilterParameters1 = subFilter1.FilterParameterList();
    EXPECT_EQ(1, subFilterParameters1.Size());

    const FilterParameter& param1 = subFilterParameters1[0];
    EXPECT_EQ(string("rq"),    param1.Name());
    EXPECT_EQ(string(">0.85"), param1.Value());

    const SubDataSet& subdataset2 = subdatasets[1];
    EXPECT_EQ(string("Alignments to chromosome 1"),           subdataset2.Name());
    EXPECT_EQ(string("ac95d0a3-94b8-4918-b3af-a3f81bbe519c"), subdataset2.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                subdataset2.Version());

    const Filters& subFilters2 = subdataset2.FilterList();
    EXPECT_EQ(1, subFilters2.Size());

    const Filter& subFilter2 = subFilters2[0];
    const FilterParameters& subFilterParameters2 = subFilter2.FilterParameterList();
    EXPECT_EQ(1, subFilterParameters2.Size());

    const FilterParameter& param2 = subFilterParameters2[0];
    EXPECT_EQ(string("RNAME"), param2.Name());
    EXPECT_EQ(string("chr1"),  param2.Value());
}

static void TestExampleXml(void)
{
    const DataSet dataset(exampleXmlFn);

    // header attributes
    EXPECT_EQ(DataSetType::SUBREADSET, dataset.Type());
    EXPECT_EQ(string(), dataset.Name());
    EXPECT_EQ(string(), dataset.Tags());
    EXPECT_EQ(string(), dataset.UniqueId());
    EXPECT_EQ(string(), dataset.Version());

    // external refs, iterator access
    const ExternalDataReferences& refs = dataset.ExternalDataReferenceList();
    EXPECT_EQ(2, dataset.NumExternalDataReferences());
    EXPECT_EQ(2, refs.Size());
    size_t i = 0;
    for (auto ref : refs) {
        if (i == 0) {
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"),                   ref.MetaType());
            EXPECT_EQ(string("file:///Users/derek/development/data/ex1.bam"), ref.ResourceId());
        }
        else if (i == 1) {
            EXPECT_EQ(string("SubreadFile.SubreadBamFile"),                   ref.MetaType());
            EXPECT_EQ(string("file:///Users/derek/development/data/ex2.bam"), ref.ResourceId());
        }
        ++i;
    }

    // filters
    const Filters& filters = dataset.FilterList();
    EXPECT_EQ(1, dataset.NumFilters());
    EXPECT_EQ(1, filters.Size());
    for ( auto filter : filters ) {

        const FilterParameters& params = filter.FilterParameterList();
        EXPECT_EQ(1, params.Size());
        for ( auto parameter : params ) {
            EXPECT_EQ(string("rq"),    parameter.Name());
            EXPECT_EQ(string(">0.75"), parameter.Value());
        }
    }

    // subdatasets
    EXPECT_EQ(1, dataset.NumSubDataSets());
    for ( auto subDataSet : dataset.SubDataSetList() ) {

        EXPECT_EQ(string("High quality"), subDataSet.Name());

        const Filters& filters = subDataSet.FilterList();
        EXPECT_EQ(1, filters.Size());
        for ( auto filter : filters ) {

            const FilterParameters& params = filter.FilterParameterList();
            EXPECT_EQ(1, params.Size());
            for ( auto parameter : params ) {
                EXPECT_EQ(string("rq"),    parameter.Name());
                EXPECT_EQ(string(">0.85"), parameter.Value());
            }
        }
    }

    // metadata
    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("5000"), metadata.TotalLength());
    EXPECT_EQ(string("500"),  metadata.NumRecords());
}

static void TestSubread1Xml(void)
{
    const DataSet dataset(subread1XmlFn);

    EXPECT_EQ(DataSetType::SUBREADSET,                        dataset.Type());
    EXPECT_EQ(string("2015-01-27T09:00:01"),                  dataset.CreatedAt());
    EXPECT_EQ(string("PacBio.DataSet.SubreadSet"),            dataset.MetaType());
    EXPECT_EQ(string("DataSet_SubreadSet"),                   dataset.Name());
    EXPECT_EQ(string("barcode moreTags mapping mytags"),      dataset.Tags());
    EXPECT_EQ(string("b095d0a3-94b8-4918-b3af-a3f81bbe519c"), dataset.UniqueId());
    EXPECT_EQ(string("2.3.0"),                                dataset.Version());

    // external refs, index access
    ExternalDataReference ref;
    const ExternalDataReferences& refs = dataset.ExternalDataReferenceList();
    EXPECT_EQ(2, refs.Size());

    ref = refs[0];
    EXPECT_EQ(string("Points to an example Subreads BAM file."), ref.Description());
    EXPECT_EQ(string("First Subreads BAM"),              ref.Name());
    EXPECT_EQ(string("SubreadFile.SubreadBamFile"),      ref.MetaType());
    EXPECT_EQ(string("file:/mnt/path/to/subreads0.bam"), ref.ResourceId());
    EXPECT_EQ(string("Example"),                         ref.Tags());

    ref = refs[1];
    EXPECT_EQ(string("Points to another example Subreads BAM file."), ref.Description());
    EXPECT_EQ(string("Second Subreads BAM"),             ref.Name());
    EXPECT_EQ(string("SubreadFile.SubreadBamFile"),      ref.MetaType());
    EXPECT_EQ(string("file:/mnt/path/to/subreads1.bam"), ref.ResourceId());
    EXPECT_EQ(string("Example"),                         ref.Tags());

    // filters
    Filter f;
    FilterParameters params;
    const Filters& filters = dataset.FilterList();
    EXPECT_EQ(2, filters.Size());

    f = filters[0];
    params = f.FilterParameterList();
    EXPECT_EQ(1, params.Size());

    EXPECT_EQ(string("rq"),    params[0].Name());
    EXPECT_EQ(string(">0.75"), params[0].Value());

    f = filters[1];
    params = f.FilterParameterList();
    EXPECT_EQ(1, params.Size());

    EXPECT_EQ(string("QNAME"),       params[0].Name());
    EXPECT_EQ(string("100/0/0_100"), params[0].Value());

    // subdatasets
    EXPECT_EQ(0, dataset.NumSubDataSets());

    // metadata
    const DataSetMetadata& metadata = dataset.Metadata();
    EXPECT_EQ(string("500000"), metadata.TotalLength());
    EXPECT_EQ(string("500"),    metadata.NumRecords());
}

TEST(DataSetIOTest, ToXml)
{
    // top-level data

    DataSet dataset(DataSetType::ALIGNMENTSET);
    dataset.CreatedAt("2015-01-27T09:00:01")
           .MetaType("PacBio.DataSet.AlignmentSet")
           .Name("DataSet_AlignmentSet")
           .Tags("barcode moreTags mapping mytags")
           .UniqueId("b095d0a3-94b8-4918-b3af-a3f81bbe519c")
           .Version("2.3.0");

    dataset.Attribute("xmlns","http://pacificbiosciences.com/PacBioDataModel.xsd");
    dataset.Attribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    dataset.Attribute("xsi:schemaLocation",
        "http://pacificbiosciences.com/PacBioDataModel.xsd "
        "../../../../../../../../common/datamodel/SequEl/EndToEnd/xsd/PacBioSecondaryDataModel.xsd");

    // external data references

    ExternalDataReference ref1;
    ref1.Name("Third Alignments BAM")
        .Description("Points to an example Alignments BAM file.")
        .MetaType("AlignmentFile.AlignmentBamFile")
        .ResourceId("file:/mnt/path/to/alignments2.bam")
        .Tags("Example");

    internal::DataSetElement pbi1("PacBioIndex");
    pbi1.Attribute("ResourceId", "file:/mnt/path/to/alignments2.pbi");
    ref1.AddChild(pbi1);

    dataset.AddExternalDataReference(ref1);

    ExternalDataReference ref2;
    ref2.Name("Fourth Alignments BAM")
        .Description("Points to another example Alignments BAM file, by relative path.")
        .MetaType("AlignmentFile.AlignmentBamFile")
        .ResourceId("file:./alignments3.bam")
        .Tags("Example");

    internal::DataSetElement pbi2("PacBioIndex");
    pbi2.Attribute("ResourceId", "file:/mnt/path/to/alignments3.pbi");
    ref2.AddChild(pbi2);

    dataset.AddExternalDataReference(ref2);

    // sub-datasets with filters

    SubDataSet subDataSet1;
    subDataSet1.Name("HighQuality Read Alignments")
               .UniqueId("ab95d0a3-94b8-4918-b3af-a3f81bbe519c")
               .Version("2.3.0");

    FilterParameter param1;
    param1.Name("rq").Value(">0.85");

    Filter filter1;
    filter1.AddParameter(param1);

    subDataSet1.AddFilter(filter1);
    dataset.AddSubDataSet(subDataSet1);

    SubDataSet subDataSet2;
    subDataSet2.Name("Alignments to chromosome 1")
               .UniqueId("ac95d0a3-94b8-4918-b3af-a3f81bbe519c")
               .Version("2.3.0");

    FilterParameter param2;
    param2.Name("RNAME").Value("chr1");

    Filter filter2;
    filter2.AddParameter(param2);

    subDataSet2.AddFilter(filter2);
    dataset.AddSubDataSet(subDataSet2);

    // write dataset
    const string expectedXml =
        "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        "<AlignmentSet CreatedAt=\"2015-01-27T09:00:01\" MetaType=\"PacBio.DataSet.AlignmentSet\" "
                "Name=\"DataSet_AlignmentSet\" Tags=\"barcode moreTags mapping mytags\" "
                "UniqueId=\"b095d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"2.3.0\" "
                "xmlns=\"http://pacificbiosciences.com/PacBioDataModel.xsd\" "
                "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
                "xsi:schemaLocation=\"http://pacificbiosciences.com/PacBioDataModel.xsd "
                "../../../../../../../../common/datamodel/SequEl/EndToEnd/xsd/PacBioSecondaryDataModel.xsd\">\n"
        "\t<ExternalDataReferences>\n"
        "\t\t<ExternalDataReference Description=\"Points to an example Alignments BAM file.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" Name=\"Third Alignments BAM\" "
                "ResourceId=\"file:/mnt/path/to/alignments2.bam\" Tags=\"Example\">\n"
        "\t\t\t<PacBioIndex ResourceId=\"file:/mnt/path/to/alignments2.pbi\" />\n"
        "\t\t</ExternalDataReference>\n"
        "\t\t<ExternalDataReference Description=\"Points to another example Alignments BAM file, by relative path.\" "
                "MetaType=\"AlignmentFile.AlignmentBamFile\" Name=\"Fourth Alignments BAM\" "
                "ResourceId=\"file:./alignments3.bam\" Tags=\"Example\">\n"
        "\t\t\t<PacBioIndex ResourceId=\"file:/mnt/path/to/alignments3.pbi\" />\n"
        "\t\t</ExternalDataReference>\n"
        "\t</ExternalDataReferences>\n"
        "\t<DataSets>\n"
        "\t\t<DataSet Name=\"HighQuality Read Alignments\" UniqueId=\"ab95d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"2.3.0\">\n"
        "\t\t\t<Filters>\n"
        "\t\t\t\t<Filter>\n"
        "\t\t\t\t\t<Parameters>\n"
        "\t\t\t\t\t\t<Parameter Name=\"rq\" Value=\">0.85\" />\n"
        "\t\t\t\t\t</Parameters>\n"
        "\t\t\t\t</Filter>\n"
        "\t\t\t</Filters>\n"
        "\t\t</DataSet>\n"
        "\t\t<DataSet Name=\"Alignments to chromosome 1\" UniqueId=\"ac95d0a3-94b8-4918-b3af-a3f81bbe519c\" Version=\"2.3.0\">\n"
        "\t\t\t<Filters>\n"
        "\t\t\t\t<Filter>\n"
        "\t\t\t\t\t<Parameters>\n"
        "\t\t\t\t\t\t<Parameter Name=\"RNAME\" Value=\"chr1\" />\n"
        "\t\t\t\t\t</Parameters>\n"
        "\t\t\t\t</Filter>\n"
        "\t\t\t</Filters>\n"
        "\t\t</DataSet>\n"
        "\t</DataSets>\n"
        "</AlignmentSet>\n";

    stringstream s;
    dataset.WriteToStream(s);
    EXPECT_EQ(expectedXml, s.str());
}

TEST(DataSetIOTest, ThrowsOnNonexistentXmlFile)
{
    EXPECT_THROW(
    {
        DataSet dataset("does/not/exist.xml");

    }, std::exception);
}

TEST(DataSetIOTest, ThrowsOnNonexistentFofnFile)
{
//    EXPECT_THROW(
//    {
//        DataSet dataset("does/not/exist.fofn");

//    }, std::exception);
}

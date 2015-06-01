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

#include <gtest/gtest.h>
#include <pbbam/dataset/DataSet.h>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace tests {

static inline
DataSet CreateDataSet(void)
{
    DataSet d;
    d.Name("foo");
    return d;
}

} // namespace tests

TEST(DataSetCoreTest, DefaultsOk)
{
    DataSet dataset;
    EXPECT_EQ(DataSetType::GENERIC, dataset.Type());
    EXPECT_TRUE(dataset.Name().empty());
    EXPECT_TRUE(dataset.Tags().empty());
    EXPECT_TRUE(dataset.UniqueId().empty());
    EXPECT_TRUE(dataset.Version().empty());

    EXPECT_EQ(0, dataset.ExternalDataReferenceList().Size());
    EXPECT_EQ(0, dataset.FilterList().Size());
    EXPECT_EQ(0, dataset.SubDataSetList().Size());
}

TEST(DataSetCoreTest, BasicGettersSettersOk)
{
    DataSet dataset;
    dataset.Name("foo");

    EXPECT_EQ(string("foo"), dataset.Name());
}

TEST(DataSetCoreTest, CopyOk)
{
    DataSet d1;
    d1.Name("foo");

    // copy ctor
    DataSet d2(d1);
    EXPECT_EQ(string("foo"), d2.Name());

    // copy assignment
    DataSet d3;
    d3 = d1;
    EXPECT_EQ(string("foo"), d3.Name());
}

TEST(DataSetCoreTest, MoveOk)
{
    DataSet d1;
    d1.Name("foo");

    // move ctor
    DataSet d2(std::move(tests::CreateDataSet()));
    EXPECT_EQ(string("foo"), d2.Name());

    // copy assignment
    DataSet d3;
    d3 = std::move(tests::CreateDataSet());
    EXPECT_EQ(string("foo"), d3.Name());
}

TEST(DataSetCoreTest, AddExternalRefs)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumExternalDataReferences());

    ExternalDataReference ref1;
    ref1.Name("file1");

    ExternalDataReference ref2;
    ref2.Name("file2");

    dataset.AddExternalDataReference(ref1);
    dataset.AddExternalDataReference(ref2);

    EXPECT_EQ(2, dataset.NumExternalDataReferences());

    DataSet dataset2;
    ExternalDataReference reusedRef;
    reusedRef.Name("file1");
    dataset2.AddExternalDataReference(reusedRef);
    reusedRef.Name("file2");
    dataset2.AddExternalDataReference(reusedRef);
    EXPECT_EQ(2, dataset2.NumExternalDataReferences());

    const ExternalDataReferences& refs = dataset2.ExternalDataReferenceList();
    EXPECT_EQ(2, refs.Size());
    EXPECT_EQ(string("file1"), refs[0].Name());
    EXPECT_EQ(string("file2"), refs[1].Name());
    string(breakl);
}

TEST(DataSetCoreTest, EditExternalRefs)
{
    DataSet dataset;

    ExternalDataReference ref;
    ref.Name("file1");
    dataset.AddExternalDataReference(ref);

    ref.Name("file2");
    dataset.AddExternalDataReference(ref);
    EXPECT_EQ(2, dataset.NumExternalDataReferences());

    dataset.ExternalDataReferenceList()[0].Name("some new name");
    EXPECT_EQ(string("some new name"), dataset.ExternalDataReferenceList()[0].Name());
}

TEST(DataSetCoreTest, AddFilters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumFilters());

    FilterParameter param1;
    param1.Name("rq");
    param1.Value(">0.85");

    Filter filter1;
    filter1.AddParameter(param1);
    EXPECT_EQ(1, filter1.NumFilterParameters());

    dataset.AddFilter(filter1);
    dataset.AddFilter(Filter());

    EXPECT_EQ(2, dataset.NumFilters());
    EXPECT_EQ(1, dataset.FilterList()[0].NumFilterParameters());
    EXPECT_EQ(0, dataset.FilterList()[1].NumFilterParameters());

    const FilterParameter& p = dataset.FilterList()[0].FilterParameterList()[0];
    EXPECT_EQ(string("rq"),    p.Name());
    EXPECT_EQ(string(">0.85"), p.Value());
}

TEST(DataSetCoreText, EditFilters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumFilters());

    FilterParameter param1;
    param1.Name("rq");
    param1.Value(">0.85");

    Filter filter1;
    filter1.AddParameter(param1);
    EXPECT_EQ(1, filter1.NumFilterParameters());

    dataset.AddFilter(filter1);
    dataset.AddFilter(Filter());

    EXPECT_EQ(2, dataset.NumFilters());
    EXPECT_EQ(1, dataset.FilterList()[0].NumFilterParameters());
    EXPECT_EQ(0, dataset.FilterList()[1].NumFilterParameters());

    FilterParameter& p = dataset.FilterList()[0].FilterParameterList()[0];
    p.Name("newParamName");
    p.Value("newParamValue");

    FilterParameter& p2 = dataset.FilterList()[0].FilterParameterList()[0];
    EXPECT_EQ(string("newParamName"),    p2.Name());
    EXPECT_EQ(string("newParamValue"), p2.Value());
}

TEST(DataSetCoreTest, AddSubDataSets)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumSubDataSets());

    SubDataSet sub1;
    sub1.Name("subset_1");

    SubDataSet sub2;
    sub2.Name("subset_2");

    dataset.AddSubDataSet(sub1);
    dataset.AddSubDataSet(sub2);
    EXPECT_EQ(2, dataset.NumSubDataSets());
}

TEST(DataSetCoreTest, EditSubDataSets)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumSubDataSets());

    SubDataSet sub1;
    sub1.Name("subset_1");
    EXPECT_EQ(string(), sub1.CreatedAt());

    SubDataSet sub2;
    sub2.Name("subset_2");

    dataset.AddSubDataSet(sub1);
    dataset.AddSubDataSet(sub2);
    EXPECT_EQ(2, dataset.NumSubDataSets());

    dataset.SubDataSetList()[0].CreatedAt("now");
    EXPECT_EQ(string("now"), dataset.SubDataSetList()[0].CreatedAt());
}

TEST(DataSetCoreTest, RemoveExternalRefs)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumExternalDataReferences());

    ExternalDataReference ref1;
    ref1.Name("file1");

    ExternalDataReference ref2;
    ref2.Name("file2");

    dataset.AddExternalDataReference(ref1);
    dataset.AddExternalDataReference(ref2);
    EXPECT_EQ(2, dataset.NumExternalDataReferences());

    dataset.RemoveExternalDataReference(ref1);
    dataset.RemoveExternalDataReference(ref2);
    EXPECT_EQ(0, dataset.NumExternalDataReferences());
}

TEST(DataSetCoreTest, RemoveFilters)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumFilters());

    Filter filter1;
    Filter filter2;
    dataset.AddFilter(filter1);
    dataset.AddFilter(filter2);
    EXPECT_EQ(2, dataset.NumFilters());

    dataset.RemoveFilter(filter1);
    dataset.RemoveFilter(filter2);
    EXPECT_EQ(0, dataset.NumFilters());
}

TEST(DataSetCoreTest, RemoveSubDataSets)
{
    DataSet dataset;
    EXPECT_EQ(0, dataset.NumSubDataSets());

    SubDataSet sub1;
    sub1.Name("subset_1");

    SubDataSet sub2;
    sub2.Name("subset_2");

    dataset.AddSubDataSet(sub1);
    dataset.AddSubDataSet(sub2);
    EXPECT_EQ(2, dataset.NumSubDataSets());

    dataset.RemoveSubDataSet(sub1);
    dataset.RemoveSubDataSet(sub2);
    EXPECT_EQ(0, dataset.NumSubDataSets());
}

TEST(DataSetCoreTest, ConversionOk)
{
    DataSet generic;
    DataSet alignment(DataSetType::ALIGNMENTSET);
    DataSet barcode(DataSetType::BARCODESET);
    DataSet ccs(DataSetType::CCSREADSET);
    DataSet contig(DataSetType::CONTIGSET);
    DataSet reference(DataSetType::REFERENCESET);
    DataSet subread(DataSetType::SUBREADSET);

    EXPECT_NO_THROW(alignment.ToAlignmentSet());
    EXPECT_NO_THROW(barcode.ToBarcodeSet());
    EXPECT_NO_THROW(ccs.ToCcsReadSet());
    EXPECT_NO_THROW(contig.ToContigSet());
    EXPECT_NO_THROW(reference.ToReferenceSet());
    EXPECT_NO_THROW(subread.ToSubreadSet());


}

// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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
#include <pbbam/SamHeader.h>
#include <pbbam/SamHeaderCodec.h>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace tests {

struct BamHdrDeleter
{
    void operator()(bam_hdr_t* hdr) {
        if (hdr)
            bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

} // namespace tests

TEST(DictionaryTest, DefaultConstruction)
{
    SamHeader::ReadGroupDictionary readGroups;

    EXPECT_TRUE(readGroups.IsEmpty());
    EXPECT_EQ(0, readGroups.Size());
    EXPECT_FALSE(readGroups.Contains("foo"));
    EXPECT_FALSE(readGroups.Contains(SamReadGroup("fake")));
}

TEST(DictionaryTest, AddUniqueAndLookup)
{
    SamReadGroup rg1("foo");
    SamReadGroup rg2("bar");
    rg2.description = "testing";

    SamHeader::ReadGroupDictionary readGroups;
    const bool add1 = readGroups.Add(rg1);
    const bool add2 = readGroups.Add(rg2);
    const bool add3 = readGroups.Add("from_string");

    EXPECT_FALSE(readGroups.IsEmpty());
    EXPECT_EQ(3, readGroups.Size());

    EXPECT_TRUE(add1);
    EXPECT_TRUE(add2);
    EXPECT_TRUE(add3);

    EXPECT_TRUE(readGroups.Contains("foo")); // added object, lookup by key
    EXPECT_TRUE(readGroups.Contains(rg1));   // added object, lookup by object
    EXPECT_TRUE(readGroups.Contains("from_string"));   // added key, lookup by key

    EXPECT_EQ(string("testing"), readGroups["bar"].description); // operator[]
}

TEST(DictionaryTest, AddListAndLookup)
{
    SamReadGroup rg1("foo");
    SamReadGroup rg2("bar");
    rg2.description = "testing";
    SamReadGroup rg3("from_string");

    vector<SamReadGroup> readGroupList;
    readGroupList.push_back(rg1);
    readGroupList.push_back(rg2);
    readGroupList.push_back(rg3);

    SamHeader::ReadGroupDictionary readGroups;
    readGroups.Add(readGroupList);

    EXPECT_FALSE(readGroups.IsEmpty());
    EXPECT_EQ(3, readGroups.Size());

    EXPECT_TRUE(readGroups.Contains("foo")); // added object, lookup by key
    EXPECT_TRUE(readGroups.Contains(rg1));   // added object, lookup by object
    EXPECT_TRUE(readGroups.Contains("from_string"));   // added key, lookup by key

    EXPECT_EQ(string("testing"), readGroups["bar"].description); // operator[]
}

TEST(DictionaryTest, AddDuplicateNotAllowed)
{
    SamReadGroup rg1("foo");
    SamReadGroup rg2("foo");
    rg2.description = "testing";

    SamHeader::ReadGroupDictionary readGroups;
    const bool add1 = readGroups.Add(rg1);
    const bool add2 = readGroups.Add(rg2);
    const bool add3 = readGroups.Add("from_string");
    const bool add4 = readGroups.Add("foo");

    EXPECT_FALSE(readGroups.IsEmpty());
    EXPECT_EQ(2, readGroups.Size());

    EXPECT_TRUE(add1);
    EXPECT_FALSE(add2);
    EXPECT_TRUE(add3);
    EXPECT_FALSE(add4);

    EXPECT_TRUE(readGroups.Contains("foo")); // added object, lookup by key
    EXPECT_TRUE(readGroups.Contains(rg1));   // added object, lookup by object
    EXPECT_TRUE(readGroups.Contains("from_string"));   // added key, lookup by key
}

TEST(DictionaryTest, IterationOk)
{
    SamReadGroup rg1("foo");
    SamReadGroup rg2("bar");
    SamReadGroup rg3("baz");

    SamHeader::ReadGroupDictionary readGroups;
    readGroups.Add(rg1);
    readGroups.Add(rg2);
    readGroups.Add(rg3);

    EXPECT_EQ(3, readGroups.Size());

    SamHeader::ReadGroupDictionary::ConstIterator iter = readGroups.ConstBegin();
    SamHeader::ReadGroupDictionary::ConstIterator end  = readGroups.ConstEnd();
    for ( ; iter != end; ++iter ) {
        const SamReadGroup& rg = (*iter);
        EXPECT_TRUE(rg.id == string("foo") ||
                    rg.id == string("bar") ||
                    rg.id == string("baz") );
    }
}

TEST(DictionaryTest, RemoveOk)
{
    SamReadGroup rg1("foo");
    SamReadGroup rg2("bar");
    SamReadGroup rg3("baz");

    SamHeader::ReadGroupDictionary readGroups;
    readGroups.Add(rg1);
    readGroups.Add(rg2);
    readGroups.Add(rg3);

    EXPECT_EQ(3, readGroups.Size());

    bool removed2 = readGroups.Remove("bar");

    EXPECT_TRUE(removed2);
    EXPECT_EQ(2, readGroups.Size());
    EXPECT_FALSE(readGroups.Contains(rg2));
    EXPECT_FALSE(readGroups.Contains("bar"));

    bool removedDummy = readGroups.Remove("__dummy__");

    EXPECT_FALSE(removedDummy);
    EXPECT_EQ(2, readGroups.Size());
}

TEST(SamHeaderTest, DefaultConstruction)
{
    SamHeader header;
    EXPECT_TRUE(header.version.empty());
    EXPECT_TRUE(header.sortOrder.empty()); // default to unknown ?
    EXPECT_TRUE(header.readGroups.IsEmpty());
    EXPECT_TRUE(header.sequences.IsEmpty());
    EXPECT_TRUE(header.programs.IsEmpty());
    EXPECT_TRUE(header.comments.empty());
}

TEST(SamHeaderCodecTest, DecodeTest)
{
    const string& text = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                         "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                         "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                         "@RG\tID:rg1\tSM:control\n"
                         "@RG\tID:rg2\tSM:condition1\n"
                         "@RG\tID:rg3\tSM:condition1\n"
                         "@PG\tID:_foo_\tPN:ide\n"
                         "@CO\tipsum and so on\n"
                         "@CO\tcitation needed\n";

    SamHeader header = SamHeaderCodec::Decode(text);

    EXPECT_EQ(string("1.1"),       header.version);
    EXPECT_EQ(string("queryname"), header.sortOrder);
    EXPECT_EQ(string("3.0b3"),     header.pacbioBamVersion);

    EXPECT_EQ(3, header.readGroups.Size());
    EXPECT_TRUE(header.readGroups.Contains("rg1"));
    EXPECT_TRUE(header.readGroups.Contains("rg2"));
    EXPECT_TRUE(header.readGroups.Contains("rg3"));
    EXPECT_EQ(string("control"),    header.readGroups["rg1"].sample);
    EXPECT_EQ(string("condition1"), header.readGroups["rg2"].sample);
    EXPECT_EQ(string("condition1"), header.readGroups["rg3"].sample);

    EXPECT_EQ(2, header.sequences.Size());
    EXPECT_TRUE(header.sequences.Contains("chr1"));
    EXPECT_TRUE(header.sequences.Contains("chr2"));
    EXPECT_EQ(string("chocobo"), header.sequences["chr1"].species);
    EXPECT_EQ(string("chocobo"), header.sequences["chr2"].species);
    EXPECT_EQ(string("2038"), header.sequences["chr1"].length);
    EXPECT_EQ(string("3042"), header.sequences["chr2"].length);

    EXPECT_EQ(1, header.programs.Size());
    EXPECT_TRUE(header.programs.Contains("_foo_"));
    EXPECT_EQ(string("ide"), header.programs["_foo_"].name);

    EXPECT_EQ(2, header.comments.size());
    EXPECT_EQ(string("ipsum and so on"), header.comments.at(0));
    EXPECT_EQ(string("citation needed"), header.comments.at(1));
}

TEST(SamHeaderCodecTest, EncodeTest)
{
    SamHeader header;
    header.version = "1.1";
    header.sortOrder = "queryname";
    header.pacbioBamVersion = "3.0b3";
    header.readGroups["rg1"].sample = "control";
    header.readGroups["rg2"].sample = "condition1";
    header.readGroups["rg3"].sample = "condition1";
    header.sequences["chr1"].length = "2038";
    header.sequences["chr1"].species = "chocobo";
    header.sequences["chr2"].length = "3042";
    header.sequences["chr2"].species = "chocobo";
    header.programs["_foo_"].name = "ide";
    header.comments.push_back("ipsum and so on");
    header.comments.push_back("citation needed");

    const string& expectedText = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                                 "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                                 "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                                 "@RG\tID:rg1\tSM:control\n"
                                 "@RG\tID:rg2\tSM:condition1\n"
                                 "@RG\tID:rg3\tSM:condition1\n"
                                 "@PG\tID:_foo_\tPN:ide\n"
                                 "@CO\tipsum and so on\n"
                                 "@CO\tcitation needed\n";

    const string& text = SamHeaderCodec::Encode(header);
    EXPECT_EQ(expectedText, text);
}

TEST(SamHeaderTest, ConvertToRawDataOk)
{
    SamHeader header;
    header.version = "1.1";
    header.sortOrder = "queryname";
    header.pacbioBamVersion = "3.0b3";
    header.readGroups["rg1"].sample = "control";
    header.readGroups["rg2"].sample = "condition1";
    header.readGroups["rg3"].sample = "condition1";
    header.sequences["chr1"].length = "2038";
    header.sequences["chr1"].species = "chocobo";
    header.sequences["chr2"].length = "3042";
    header.sequences["chr2"].species = "chocobo";
    header.programs["_foo_"].name = "ide";
    header.comments.push_back("ipsum and so on");
    header.comments.push_back("citation needed");

    const string& expectedText = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                                 "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                                 "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                                 "@RG\tID:rg1\tSM:control\n"
                                 "@RG\tID:rg2\tSM:condition1\n"
                                 "@RG\tID:rg3\tSM:condition1\n"
                                 "@PG\tID:_foo_\tPN:ide\n"
                                 "@CO\tipsum and so on\n"
                                 "@CO\tcitation needed\n";

    shared_ptr<bam_hdr_t> rawData = header.CreateRawData();
    const string rawText = string(rawData->text, rawData->l_text);
    EXPECT_EQ(expectedText, rawText);
}

TEST(SamHeaderTest, ExtractFromRawDataOk)
{
    SamHeader header;
    header.version = "1.1";
    header.sortOrder = "queryname";
    header.pacbioBamVersion = "3.0b3";
    header.readGroups["rg1"].sample = "control";
    header.readGroups["rg2"].sample = "condition1";
    header.readGroups["rg3"].sample = "condition1";
    header.sequences["chr1"].length = "2038";
    header.sequences["chr1"].species = "chocobo";
    header.sequences["chr2"].length = "3042";
    header.sequences["chr2"].species = "chocobo";
    header.programs["_foo_"].name = "ide";
    header.comments.push_back("ipsum and so on");
    header.comments.push_back("citation needed");

    const string& expectedText = "@HD\tVN:1.1\tSO:queryname\tpb:3.0b3\n"
                                 "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                                 "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                                 "@RG\tID:rg1\tSM:control\n"
                                 "@RG\tID:rg2\tSM:condition1\n"
                                 "@RG\tID:rg3\tSM:condition1\n"
                                 "@PG\tID:_foo_\tPN:ide\n"
                                 "@CO\tipsum and so on\n"
                                 "@CO\tcitation needed\n";

    shared_ptr<bam_hdr_t> rawData = header.CreateRawData();
    SamHeader newHeader = SamHeader::FromRawData(rawData);
    EXPECT_EQ(header.version, newHeader.version);
    EXPECT_EQ(header.sortOrder, newHeader.sortOrder);
    EXPECT_EQ(header.pacbioBamVersion, newHeader.pacbioBamVersion);

    const string& text = SamHeaderCodec::Encode(newHeader);
    EXPECT_EQ(expectedText, text);
}

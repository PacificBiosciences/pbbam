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

#include <gtest/gtest.h>

#ifdef PBBAM_TESTING
#define private public
#endif

#include <pbbam/BamRecordImpl.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

// NOTE: these tests check "high-level" tag query/manipulation via BamRecordImpl.
//       For raw Tag/TagCollection tests, see test_Tags.cpp
//       For encoding tests, see test_BamRecordImplVariableData.cpp

TEST(BamRecordImplTagsTest, HasTagTest)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    EXPECT_FALSE(bam.HasTag("zz"));
    EXPECT_FALSE(bam.HasTag(""));
    EXPECT_FALSE(bam.HasTag("some_too_long_name"));

    const TagCollection& fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_TRUE(fetchedTags.Contains("XY"));
    EXPECT_FALSE(fetchedTags.Contains("zz"));
    EXPECT_FALSE(fetchedTags.Contains(""));
    EXPECT_FALSE(fetchedTags.Contains("some_too_long_name"));
}

TEST(BamRecordImplTagsTest, SimpleAddTag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_FALSE(bam.HasTag("XY"));

    bam.AddTag("XY", (int32_t)-42);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    const TagCollection& fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_TRUE(fetchedTags.Contains("XY"));
    EXPECT_FALSE(fetchedTags.Contains("zz"));
    EXPECT_FALSE(fetchedTags.Contains(""));
    EXPECT_FALSE(fetchedTags.Contains("some_too_long_name"));

    EXPECT_EQ(-42, fetchedTags.at("XY").ToInt32());

    // fail on invalid adds
    EXPECT_FALSE(bam.AddTag("", (int32_t)-42));
    EXPECT_FALSE(bam.AddTag("some_too_long_name", (int32_t)-42));
    EXPECT_FALSE(bam.AddTag("XY", (int32_t)-42));                   // reject duplicate
}

TEST(BamRecordImplTagsTest, SimpleRemoveTag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    const bool removedOk = bam.RemoveTag("XY");
    EXPECT_TRUE(removedOk);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_FALSE(bam.HasTag("XY"));

    const TagCollection& fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_FALSE(fetchedTags.Contains("XY"));
    EXPECT_FALSE(fetchedTags.Contains("zz"));
    EXPECT_FALSE(fetchedTags.Contains(""));
    EXPECT_FALSE(fetchedTags.Contains("some_too_long_name"));

    // fail on invalid removes
    EXPECT_FALSE(bam.RemoveTag(""));
    EXPECT_FALSE(bam.RemoveTag("some_too_long_name"));
    EXPECT_FALSE(bam.RemoveTag("zz"));                 // reject remove unknown
}

TEST(BamRecordImplTagsTest, SimpleEditTag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("XY"));

    const TagCollection& fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_TRUE(fetchedTags.Contains("XY"));
    EXPECT_EQ(-42, fetchedTags.at("XY").ToInt32());

    const bool editedOk = bam.EditTag("XY", (int32_t)500);
    EXPECT_TRUE(editedOk);
    EXPECT_TRUE(bam.HasTag("XY"));

    const TagCollection& fetchedTags2 = bam.Tags();
    EXPECT_TRUE(fetchedTags2.Contains("HX"));
    EXPECT_TRUE(fetchedTags2.Contains("CA"));
    EXPECT_TRUE(fetchedTags2.Contains("XY"));
    EXPECT_EQ(500, fetchedTags2.at("XY").ToInt32());

    // fail on invalid edits
    EXPECT_FALSE(bam.EditTag("", 500));
    EXPECT_FALSE(bam.EditTag("some_too_long_name", 500));
    EXPECT_FALSE(bam.EditTag("zz", 500));                 // reject edit unknown
}

TEST(BamRecordImplTagsTest, SimpleQueryTag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = (int32_t)-42;

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("XY"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    EXPECT_EQ(string("1abc75"),              bam.TagValue("HX").ToString());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), bam.TagValue("CA").ToUInt8Array());
    EXPECT_EQ((int32_t)-42,                  bam.TagValue("XY").ToInt32());

    EXPECT_FALSE(bam.HasTag("zz"));
    EXPECT_FALSE(bam.HasTag(""));
    EXPECT_FALSE(bam.HasTag("some_too_long_name"));

    EXPECT_EQ(Tag(), bam.TagValue("zz"));
    EXPECT_EQ(Tag(), bam.TagValue(""));
    EXPECT_EQ(Tag(), bam.TagValue("some_too_long_name"));
}


#include <pbbam/BamRecordImpl.h>

#include <cstdint>

#include <string>
#include <vector>

#include <gtest/gtest.h>

using namespace PacBio;
using namespace PacBio::BAM;

// NOTE: these tests check "high-level" tag query/manipulation via BamRecordImpl.
//       For raw Tag/TagCollection tests, see test_Tags.cpp
//       For encoding tests, see test_BamRecordImplVariableData.cpp

TEST(BAM_BamRecordImplTags, can_check_presence_of_tag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<std::uint8_t>({34, 5, 125});
    tags["XY"] = std::int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    EXPECT_FALSE(bam.HasTag("zz"));
    EXPECT_FALSE(bam.HasTag(""));
    EXPECT_FALSE(bam.HasTag("some_too_long_name"));

    const TagCollection fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_TRUE(fetchedTags.Contains("XY"));
    EXPECT_FALSE(fetchedTags.Contains("zz"));
    EXPECT_FALSE(fetchedTags.Contains(""));
    EXPECT_FALSE(fetchedTags.Contains("some_too_long_name"));
}

TEST(BAM_BamRecordImplTags, can_add_tag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<std::uint8_t>({34, 5, 125});

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_FALSE(bam.HasTag("XY"));

    bam.AddTag("XY", std::int32_t{-42});

    EXPECT_TRUE(bam.HasTag("HX"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    const TagCollection fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_TRUE(fetchedTags.Contains("XY"));
    EXPECT_FALSE(fetchedTags.Contains("zz"));
    EXPECT_FALSE(fetchedTags.Contains(""));
    EXPECT_FALSE(fetchedTags.Contains("some_too_long_name"));

    EXPECT_EQ(-42, fetchedTags.at("XY").ToInt32());

    // fail on invalid adds
    EXPECT_FALSE(bam.AddTag("", std::int32_t{-42}));
    EXPECT_FALSE(bam.AddTag("some_too_long_name", std::int32_t{-42}));
    EXPECT_FALSE(bam.AddTag("XY", std::int32_t{-42}));  // reject duplicate
}

TEST(BAM_BamRecordImplTags, can_remove_tag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<std::uint8_t>({34, 5, 125});
    tags["XY"] = std::int32_t{-42};

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

    const TagCollection fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_FALSE(fetchedTags.Contains("XY"));
    EXPECT_FALSE(fetchedTags.Contains("zz"));
    EXPECT_FALSE(fetchedTags.Contains(""));
    EXPECT_FALSE(fetchedTags.Contains("some_too_long_name"));

    // fail on invalid removes
    EXPECT_FALSE(bam.RemoveTag(""));
    EXPECT_FALSE(bam.RemoveTag("some_too_long_name"));
    EXPECT_FALSE(bam.RemoveTag("zz"));  // reject remove unknown
}

TEST(BAM_BamRecordImplTags, can_edit_tag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<std::uint8_t>({34, 5, 125});
    tags["XY"] = std::int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("XY"));

    const TagCollection fetchedTags = bam.Tags();
    EXPECT_TRUE(fetchedTags.Contains("HX"));
    EXPECT_TRUE(fetchedTags.Contains("CA"));
    EXPECT_TRUE(fetchedTags.Contains("XY"));
    EXPECT_EQ(-42, fetchedTags.at("XY").ToInt32());

    const bool editedOk = bam.EditTag("XY", std::int32_t{500});
    EXPECT_TRUE(editedOk);
    EXPECT_TRUE(bam.HasTag("XY"));

    const TagCollection fetchedTags2 = bam.Tags();
    EXPECT_TRUE(fetchedTags2.Contains("HX"));
    EXPECT_TRUE(fetchedTags2.Contains("CA"));
    EXPECT_TRUE(fetchedTags2.Contains("XY"));
    EXPECT_EQ(500, fetchedTags2.at("XY").ToInt32());

    // fail on invalid edits
    EXPECT_FALSE(bam.EditTag("", 500));
    EXPECT_FALSE(bam.EditTag("some_too_long_name", 500));
    EXPECT_FALSE(bam.EditTag("zz", 500));  // reject edit unknown
}

TEST(BAM_BamRecordImplTags, can_query_tag)
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<std::uint8_t>({34, 5, 125});
    tags["XY"] = std::int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);

    EXPECT_TRUE(bam.HasTag("XY"));
    EXPECT_TRUE(bam.HasTag("CA"));
    EXPECT_TRUE(bam.HasTag("XY"));

    EXPECT_EQ(std::string("1abc75"), bam.TagValue("HX").ToString());
    EXPECT_EQ(std::vector<std::uint8_t>({34, 5, 125}), bam.TagValue("CA").ToUInt8Array());
    EXPECT_EQ(std::int32_t{-42}, bam.TagValue("XY").ToInt32());

    EXPECT_FALSE(bam.HasTag("zz"));
    EXPECT_FALSE(bam.HasTag(""));
    EXPECT_FALSE(bam.HasTag("some_too_long_name"));

    EXPECT_EQ(Tag(), bam.TagValue("zz"));
    EXPECT_EQ(Tag(), bam.TagValue(""));
    EXPECT_EQ(Tag(), bam.TagValue("some_too_long_name"));
}

TEST(BAM_BamRecordImplTags, can_query_array_or_string_tag_length)
{
    TagCollection tags;
    tags["xx"] = std::vector<std::int32_t>{};
    tags["yy"] = std::string{};
    tags["aa"] = std::vector<std::uint8_t>({34, 5, 125});
    tags["bb"] = std::string{"triforce"};
    tags["cc"] = std::int32_t{-42};

    BamRecordImpl bam;
    bam.Tags(tags);

    // empty array
    ASSERT_TRUE(bam.TagLength("xx").has_value());
    EXPECT_EQ(0, bam.TagLength("xx").value());

    // empty string
    ASSERT_TRUE(bam.TagLength("yy").has_value());
    EXPECT_EQ(0, bam.TagLength("yy").value());

    // non-empty array
    ASSERT_TRUE(bam.TagLength("aa").has_value());
    EXPECT_EQ(3, bam.TagLength("aa").value());

    // non-empty string
    ASSERT_TRUE(bam.TagLength("bb").has_value());
    EXPECT_EQ(8, bam.TagLength("bb").value());

    // scalar value
    EXPECT_FALSE(bam.TagLength("cc").has_value());

    // does not exist
    EXPECT_FALSE(bam.TagLength("dd").has_value());
}

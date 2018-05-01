// Author: Derek Barnett

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <typeinfo>
#include <vector>

#include <gtest/gtest.h>
#include <boost/type_traits/is_convertible.hpp>

#include <pbbam/BamTagCodec.h>
#include <pbbam/SamTagCodec.h>
#include <pbbam/TagCollection.h>

using namespace PacBio;
using namespace PacBio::BAM;

TEST(TagTest, TagConstruction)
{
    int8_t i8 = 0;
    uint8_t u8 = 0;
    int16_t i16 = 0;
    uint16_t u16 = 0;
    int32_t i32 = 0;
    uint32_t u32 = 0;
    float f = 0.0;
    std::string str = "";
    std::vector<int8_t> i8_array;
    std::vector<uint8_t> u8_array;
    std::vector<int16_t> i16_array;
    std::vector<uint16_t> u16_array;
    std::vector<int32_t> i32_array;
    std::vector<uint32_t> u32_Array;
    std::vector<float> float_array;

    signed char c = 'A';
    unsigned char uc = 'A';

    Tag i8Tag(i8);
    Tag u8Tag(u8);
    Tag i16Tag(i16);
    Tag u16Tag(u16);
    Tag i32Tag(i32);
    Tag u32Tag(u32);
    Tag floatTag(f);
    Tag stringTag(str);
    Tag i8_array_Tag(i8_array);
    Tag u8_array_Tag(u8_array);
    Tag i16_array_Tag(i16_array);
    Tag u16_array_Tag(u16_array);
    Tag i32_array_Tag(i32_array);
    Tag u32_array_Tag(u32_Array);
    Tag float_array_Tag(float_array);

    Tag charTag(c, TagModifier::ASCII_CHAR);
    Tag ucharTag(uc, TagModifier::ASCII_CHAR);

    EXPECT_TRUE(i8Tag.Type() == TagDataType::INT8);
    EXPECT_TRUE(u8Tag.Type() == TagDataType::UINT8);
    EXPECT_TRUE(i16Tag.Type() == TagDataType::INT16);
    EXPECT_TRUE(u16Tag.Type() == TagDataType::UINT16);
    EXPECT_TRUE(i32Tag.Type() == TagDataType::INT32);
    EXPECT_TRUE(u32Tag.Type() == TagDataType::UINT32);
    EXPECT_TRUE(floatTag.Type() == TagDataType::FLOAT);
    EXPECT_TRUE(stringTag.Type() == TagDataType::STRING);
    EXPECT_TRUE(i8_array_Tag.Type() == TagDataType::INT8_ARRAY);
    EXPECT_TRUE(u8_array_Tag.Type() == TagDataType::UINT8_ARRAY);
    EXPECT_TRUE(i16_array_Tag.Type() == TagDataType::INT16_ARRAY);
    EXPECT_TRUE(u16_array_Tag.Type() == TagDataType::UINT16_ARRAY);
    EXPECT_TRUE(i32_array_Tag.Type() == TagDataType::INT32_ARRAY);
    EXPECT_TRUE(u32_array_Tag.Type() == TagDataType::UINT32_ARRAY);
    EXPECT_TRUE(float_array_Tag.Type() == TagDataType::FLOAT_ARRAY);

    EXPECT_TRUE(charTag.ToAscii() == 'A');
    EXPECT_TRUE(ucharTag.ToAscii() == 'A');
}

TEST(TagTest, CopyAndCompare)
{
    int8_t i8 = 0;
    uint8_t u8 = 0;
    int16_t i16 = 0;
    uint16_t u16 = 0;
    int32_t i32 = 0;
    uint32_t u32 = 0;
    float f = 0.0;
    std::string str = "";
    std::vector<int8_t> i8_array;
    std::vector<uint8_t> u8_array;
    std::vector<int16_t> i16_array;
    std::vector<uint16_t> u16_array;
    std::vector<int32_t> i32_array;
    std::vector<uint32_t> u32_Array;
    std::vector<float> float_array;

    Tag i8Tag(i8);
    Tag u8Tag(u8);
    Tag i16Tag(i16);
    Tag u16Tag(u16);
    Tag i32Tag(i32);
    Tag u32Tag(u32);
    Tag floatTag(f);
    Tag stringTag(str);
    Tag i8_array_Tag(i8_array);
    Tag u8_array_Tag(u8_array);
    Tag i16_array_Tag(i16_array);
    Tag u16_array_Tag(u16_array);
    Tag i32_array_Tag(i32_array);
    Tag u32_array_Tag(u32_Array);
    Tag float_array_Tag(float_array);

    Tag i8Tag2 = i8Tag;
    Tag u8Tag2 = u8Tag;
    Tag i16Tag2 = i16Tag;
    Tag u16Tag2 = u16Tag;
    Tag i32Tag2 = i32Tag;
    Tag u32Tag2 = u32Tag;
    Tag floatTag2 = floatTag;
    Tag stringTag2 = stringTag;
    Tag i8_array_Tag2 = i8_array_Tag;
    Tag u8_array_Tag2 = u8_array_Tag;
    Tag i16_array_Tag2 = i16_array_Tag;
    Tag u16_array_Tag2 = u16_array_Tag;
    Tag i32_array_Tag2 = i32_array_Tag;
    Tag u32_array_Tag2 = u32_array_Tag;
    Tag float_array_Tag2 = float_array_Tag;

    EXPECT_EQ(i8Tag, i8Tag2);
    EXPECT_EQ(u8Tag, u8Tag2);
    EXPECT_EQ(i16Tag, i16Tag2);
    EXPECT_EQ(u16Tag, u16Tag2);
    EXPECT_EQ(i32Tag, i32Tag2);
    EXPECT_EQ(u32Tag, u32Tag2);
    EXPECT_EQ(floatTag, floatTag2);
    EXPECT_EQ(stringTag, stringTag2);
    EXPECT_EQ(i8_array_Tag, i8_array_Tag2);
    EXPECT_EQ(u8_array_Tag, u8_array_Tag2);
    EXPECT_EQ(i16_array_Tag, i16_array_Tag2);
    EXPECT_EQ(u16_array_Tag, u16_array_Tag2);
    EXPECT_EQ(i32_array_Tag, i32_array_Tag2);
    EXPECT_EQ(u32_array_Tag, u32_array_Tag2);
    EXPECT_EQ(float_array_Tag, float_array_Tag2);
}

TEST(TagTest, Type_None)
{
    Tag tag;

    EXPECT_TRUE(tag.Type() == TagDataType::INVALID);
    EXPECT_TRUE(tag.IsNull());
    EXPECT_TRUE(tag.Typename() == "none");

    EXPECT_FALSE(tag.IsNumeric());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());
}

TEST(TagTest, Type_Int8)
{
    const int8_t v = -42;
    const Tag tag(v);

    int8_t v2{};
    EXPECT_NO_THROW(v2 = tag.ToInt8());

    EXPECT_TRUE(tag.Type() == TagDataType::INT8);
    EXPECT_TRUE(tag.Typename() == "int8_t");
    EXPECT_TRUE(tag.IsInt8());

    EXPECT_TRUE(tag.IsSignedInt());
    EXPECT_TRUE(tag.IsIntegral());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsUnsignedInt());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_UInt8)
{
    const uint8_t v = 42;
    const Tag tag(v);

    uint8_t v2{};
    EXPECT_NO_THROW(v2 = tag.ToUInt8());

    EXPECT_TRUE(tag.Type() == TagDataType::UINT8);
    EXPECT_TRUE(tag.Typename() == "uint8_t");
    EXPECT_TRUE(tag.IsUInt8());

    EXPECT_TRUE(tag.IsUnsignedInt());
    EXPECT_TRUE(tag.IsIntegral());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsSignedInt());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_Ascii)
{
    const char c = '$';
    const signed char sc = '$';
    const unsigned char uc = '$';
    const uint8_t u8 = 65;
    const int8_t i8 = 66;

    {  // old style: construct-then-modify

        Tag fromPlainChar = Tag(c);
        Tag fromSignedChar = Tag(sc);
        Tag fromUnsignedChar = Tag(uc);
        Tag fromUint8 = Tag(u8);
        Tag fromInt8 = Tag(i8);
        fromPlainChar.Modifier(TagModifier::ASCII_CHAR);
        fromSignedChar.Modifier(TagModifier::ASCII_CHAR);
        fromUnsignedChar.Modifier(TagModifier::ASCII_CHAR);
        fromUint8.Modifier(TagModifier::ASCII_CHAR);
        fromInt8.Modifier(TagModifier::ASCII_CHAR);

        EXPECT_TRUE(fromPlainChar.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromPlainChar.IsIntegral());
        EXPECT_TRUE(fromPlainChar.IsNumeric());
        EXPECT_EQ('$', fromPlainChar.ToAscii());

        EXPECT_TRUE(fromSignedChar.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromSignedChar.IsIntegral());
        EXPECT_TRUE(fromSignedChar.IsNumeric());
        EXPECT_EQ('$', fromSignedChar.ToAscii());

        EXPECT_TRUE(fromUnsignedChar.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromUnsignedChar.IsIntegral());
        EXPECT_TRUE(fromUnsignedChar.IsNumeric());
        EXPECT_EQ('$', fromUnsignedChar.ToAscii());

        EXPECT_TRUE(fromUint8.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromUint8.IsIntegral());
        EXPECT_TRUE(fromUint8.IsNumeric());
        EXPECT_EQ('A', fromUint8.ToAscii());

        EXPECT_TRUE(fromInt8.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromInt8.IsIntegral());
        EXPECT_TRUE(fromInt8.IsNumeric());
        EXPECT_EQ('B', fromInt8.ToAscii());
    }

    {  // new style: construct directly as ASCII

        const Tag fromPlainChar = Tag(c, TagModifier::ASCII_CHAR);
        const Tag fromSignedChar = Tag(sc, TagModifier::ASCII_CHAR);
        const Tag fromUnsignedChar = Tag(uc, TagModifier::ASCII_CHAR);
        const Tag fromUint8 = Tag(u8, TagModifier::ASCII_CHAR);
        const Tag fromInt8 = Tag(i8, TagModifier::ASCII_CHAR);

        EXPECT_TRUE(fromPlainChar.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromPlainChar.IsIntegral());
        EXPECT_TRUE(fromPlainChar.IsNumeric());
        EXPECT_EQ('$', fromPlainChar.ToAscii());

        EXPECT_TRUE(fromSignedChar.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromSignedChar.IsIntegral());
        EXPECT_TRUE(fromSignedChar.IsNumeric());
        EXPECT_EQ('$', fromSignedChar.ToAscii());

        EXPECT_TRUE(fromUnsignedChar.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromUnsignedChar.IsIntegral());
        EXPECT_TRUE(fromUnsignedChar.IsNumeric());
        EXPECT_EQ('$', fromUnsignedChar.ToAscii());

        EXPECT_TRUE(fromUint8.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromUint8.IsIntegral());
        EXPECT_TRUE(fromUint8.IsNumeric());
        EXPECT_EQ('A', fromUint8.ToAscii());

        EXPECT_TRUE(fromInt8.HasModifier(TagModifier::ASCII_CHAR));
        EXPECT_TRUE(fromInt8.IsIntegral());
        EXPECT_TRUE(fromInt8.IsNumeric());
        EXPECT_EQ('B', fromInt8.ToAscii());
    }

    // check invalid constructs
    EXPECT_THROW(Tag('A', TagModifier::HEX_STRING), std::runtime_error);
}

TEST(TagTest, Type_Int16)
{
    const int16_t v = -42;
    const Tag tag(v);

    int16_t v2{};
    EXPECT_NO_THROW(v2 = tag.ToInt16());

    EXPECT_TRUE(tag.Type() == TagDataType::INT16);
    EXPECT_TRUE(tag.Typename() == "int16_t");
    EXPECT_TRUE(tag.IsInt16());
    EXPECT_TRUE(tag.IsSignedInt());
    EXPECT_TRUE(tag.IsIntegral());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsUnsignedInt());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_UInt16)
{
    const uint16_t v = 42;
    const Tag tag(v);

    uint16_t v2;
    EXPECT_NO_THROW(v2 = tag.ToUInt16());

    EXPECT_TRUE(tag.Type() == TagDataType::UINT16);
    EXPECT_TRUE(tag.Typename() == "uint16_t");
    EXPECT_TRUE(tag.IsUInt16());
    EXPECT_TRUE(tag.IsUnsignedInt());
    EXPECT_TRUE(tag.IsIntegral());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsSignedInt());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_Int32)
{
    const int32_t v = -42;
    const Tag tag(v);

    int32_t v2;
    EXPECT_NO_THROW(v2 = tag.ToInt32());

    EXPECT_TRUE(tag.Type() == TagDataType::INT32);
    EXPECT_TRUE(tag.Typename() == "int32_t");
    EXPECT_TRUE(tag.IsInt32());
    EXPECT_TRUE(tag.IsSignedInt());
    EXPECT_TRUE(tag.IsIntegral());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsUnsignedInt());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_UInt32)
{
    const uint32_t v = 42;
    const Tag tag(v);

    uint32_t v2;
    EXPECT_NO_THROW(v2 = tag.ToUInt32());

    EXPECT_TRUE(tag.Type() == TagDataType::UINT32);
    EXPECT_TRUE(tag.Typename() == "uint32_t");
    EXPECT_TRUE(tag.IsUInt32());
    EXPECT_TRUE(tag.IsUnsignedInt());
    EXPECT_TRUE(tag.IsIntegral());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsSignedInt());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_Float)
{
    const float v = 3.141;
    const Tag tag(v);

    float v2;
    EXPECT_NO_THROW(v2 = tag.ToFloat());

    EXPECT_TRUE(tag.Type() == TagDataType::FLOAT);
    EXPECT_TRUE(tag.Typename() == "float");
    EXPECT_TRUE(tag.IsFloat());
    EXPECT_TRUE(tag.IsNumeric());

    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsIntegral());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_String)
{
    const std::string v = "foo_who";
    const Tag tag(v);

    std::string v2;
    EXPECT_NO_THROW(v2 = tag.ToString());

    EXPECT_TRUE(tag.Type() == TagDataType::STRING);
    EXPECT_TRUE(tag.Typename() == "string");
    EXPECT_TRUE(tag.IsString());

    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());
    EXPECT_FALSE(tag.IsArray());

    EXPECT_EQ(v, v2);

    // "Hex format" string
    const Tag hex("DEADBEEF", TagModifier::HEX_STRING);
    EXPECT_TRUE(hex.Type() == TagDataType::STRING);
    EXPECT_TRUE(hex.Typename() == "string");
    EXPECT_TRUE(hex.IsString());
    EXPECT_TRUE(hex.HasModifier(TagModifier::HEX_STRING));
    EXPECT_FALSE(hex.IsNull());
    EXPECT_FALSE(hex.IsNumeric());
    EXPECT_FALSE(hex.IsArray());

    // check invalid constructs
    EXPECT_THROW(Tag("DEADBEEF", TagModifier::ASCII_CHAR), std::runtime_error);
}

TEST(TagTest, Type_Int8Array)
{
    const std::vector<int8_t> v = {-42, 100, 0};
    const Tag tag(v);

    std::vector<int8_t> v2;
    EXPECT_NO_THROW(v2 = tag.ToInt8Array());

    EXPECT_TRUE(tag.Type() == TagDataType::INT8_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<int8_t>");
    EXPECT_TRUE(tag.IsInt8Array());
    EXPECT_TRUE(tag.IsSignedArray());
    EXPECT_TRUE(tag.IsIntegralArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_UInt8Array)
{
    const std::vector<uint8_t> v = {42, 200, 0};
    const Tag tag(v);

    std::vector<uint8_t> v2;
    EXPECT_NO_THROW(v2 = tag.ToUInt8Array());

    EXPECT_TRUE(tag.Type() == TagDataType::UINT8_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<uint8_t>");
    EXPECT_TRUE(tag.IsUInt8Array());
    EXPECT_TRUE(tag.IsUnsignedArray());
    EXPECT_TRUE(tag.IsIntegralArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_Int16Array)
{
    const std::vector<int16_t> v = {42, -300, 0};
    const Tag tag(v);

    std::vector<int16_t> v2;
    EXPECT_NO_THROW(v2 = tag.ToInt16Array());

    EXPECT_TRUE(tag.Type() == TagDataType::INT16_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<int16_t>");
    EXPECT_TRUE(tag.IsInt16Array());
    EXPECT_TRUE(tag.IsSignedArray());
    EXPECT_TRUE(tag.IsIntegralArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_UInt16Array)
{
    const std::vector<uint16_t> v = {42, 300, 0};
    const Tag tag(v);

    std::vector<uint16_t> v2;
    EXPECT_NO_THROW(v2 = tag.ToUInt16Array());

    EXPECT_TRUE(tag.Type() == TagDataType::UINT16_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<uint16_t>");
    EXPECT_TRUE(tag.IsUInt16Array());
    EXPECT_TRUE(tag.IsUnsignedArray());
    EXPECT_TRUE(tag.IsIntegralArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
    ;
}

TEST(TagTest, Type_Int32Array)
{
    const std::vector<int32_t> v = {42, -300, 0};
    const Tag tag(v);

    std::vector<int32_t> v2;
    EXPECT_NO_THROW(v2 = tag.ToInt32Array());

    EXPECT_TRUE(tag.Type() == TagDataType::INT32_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<int32_t>");
    EXPECT_TRUE(tag.IsInt32Array());
    EXPECT_TRUE(tag.IsSignedArray());
    EXPECT_TRUE(tag.IsIntegralArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_UInt32Array)
{
    const std::vector<uint32_t> v = {42, 300, 0};
    const Tag tag(v);

    std::vector<uint32_t> v2;
    EXPECT_NO_THROW(v2 = tag.ToUInt32Array());

    EXPECT_TRUE(tag.Type() == TagDataType::UINT32_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<uint32_t>");
    EXPECT_TRUE(tag.IsUInt32Array());
    EXPECT_TRUE(tag.IsUnsignedArray());
    EXPECT_TRUE(tag.IsIntegralArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, Type_FloatArray)
{
    const std::vector<float> v = {1.1f, 1.2f, 1.3f};
    const Tag tag(v);

    std::vector<float> v2;
    EXPECT_NO_THROW(v2 = tag.ToFloatArray());

    EXPECT_TRUE(tag.Type() == TagDataType::FLOAT_ARRAY);
    EXPECT_TRUE(tag.Typename() == "vector<float>");
    EXPECT_TRUE(tag.IsFloatArray());
    EXPECT_TRUE(tag.IsArray());

    EXPECT_FALSE(tag.IsIntegralArray());
    EXPECT_FALSE(tag.IsFloat());
    EXPECT_FALSE(tag.IsString());
    EXPECT_FALSE(tag.IsNull());
    EXPECT_FALSE(tag.IsNumeric());

    EXPECT_EQ(v, v2);
}

TEST(TagTest, CastBackToOriginalOk)
{
    int8_t i8 = 0;
    uint8_t u8 = 0;
    int16_t i16 = 0;
    uint16_t u16 = 0;
    int32_t i32 = 0;
    uint32_t u32 = 0;
    float f = 0.0;
    std::string str = "";
    std::vector<int8_t> i8_array;
    std::vector<uint8_t> u8_array;
    std::vector<int16_t> i16_array;
    std::vector<uint16_t> u16_array;
    std::vector<int32_t> i32_array;
    std::vector<uint32_t> u32_array;
    std::vector<float> float_array;

    Tag i8Tag(i8);
    Tag u8Tag(u8);
    Tag i16Tag(i16);
    Tag u16Tag(u16);
    Tag i32Tag(i32);
    Tag u32Tag(u32);
    Tag floatTag(f);
    Tag stringTag(str);
    Tag i8_array_Tag(i8_array);
    Tag u8_array_Tag(u8_array);
    Tag i16_array_Tag(i16_array);
    Tag u16_array_Tag(u16_array);
    Tag i32_array_Tag(i32_array);
    Tag u32_array_Tag(u32_array);
    Tag float_array_Tag(float_array);

    EXPECT_NO_THROW({
        i8 = i8Tag.ToInt8();
        u8 = u8Tag.ToUInt8();
        i16 = i16Tag.ToInt16();
        u16 = u16Tag.ToUInt16();
        i32 = i32Tag.ToInt32();
        u32 = u32Tag.ToUInt32();
        f = floatTag.ToFloat();
        str = stringTag.ToString();
        i8_array = i8_array_Tag.ToInt8Array();
        u8_array = u8_array_Tag.ToUInt8Array();
        i16_array = i16_array_Tag.ToInt16Array();
        u16_array = u16_array_Tag.ToUInt16Array();
        i32_array = i32_array_Tag.ToInt32Array();
        u32_array = u32_array_Tag.ToUInt32Array();
        float_array = float_array_Tag.ToFloatArray();
    });
}

TEST(TagTest, ConvertToInt8)
{
    Tag zero(int32_t{0});
    Tag min(int32_t{std::numeric_limits<int8_t>::min()});
    Tag normal(int32_t{42});
    Tag max(int32_t{std::numeric_limits<int8_t>::max()});
    Tag floatTag(float{3.14});
    Tag stringTag(std::string{"foo"});
    Tag arrayTag(std::vector<int8_t>{{1, 2, 3}});

    // allowed
    EXPECT_NO_THROW({
        zero.ToInt8();
        min.ToInt8();
        normal.ToInt8();
        max.ToInt8();
    });

    // not allowed
    EXPECT_THROW(floatTag.ToInt8(), std::exception);
    EXPECT_THROW(stringTag.ToInt8(), std::exception);
    EXPECT_THROW(arrayTag.ToInt8(), std::exception);
}

TEST(TagTest, ConvertToUInt8)
{
    Tag zero(int32_t{0});
    Tag neg(int32_t{-1});
    Tag normal(int32_t{42});
    Tag max(int32_t{std::numeric_limits<uint8_t>::max()});
    Tag floatTag(float{3.14});
    Tag stringTag(std::string{"foo"});
    Tag arrayTag(std::vector<uint8_t>{{1, 2, 3}});

    // allowed
    EXPECT_NO_THROW({
        zero.ToUInt8();
        normal.ToUInt8();
        max.ToUInt8();
    });

    // not allowed
    EXPECT_THROW(neg.ToUInt8(), std::exception);
    EXPECT_THROW(floatTag.ToUInt8(), std::exception);
    EXPECT_THROW(stringTag.ToUInt8(), std::exception);
    EXPECT_THROW(arrayTag.ToUInt8(), std::exception);
}

TEST(TagTest, ConvertToInt16)
{
    Tag zero(int32_t{0});
    Tag min(int32_t{std::numeric_limits<int16_t>::min()});
    Tag normal(int32_t{42});
    Tag max(int32_t{std::numeric_limits<int16_t>::max()});
    Tag floatTag(float{3.14});
    Tag stringTag(std::string{"foo"});
    Tag arrayTag(std::vector<int16_t>{{1, 2, 3}});

    // allowed
    EXPECT_NO_THROW({
        zero.ToInt16();
        min.ToInt16();
        normal.ToInt16();
        max.ToInt16();
    });

    // not allowed
    EXPECT_THROW(floatTag.ToInt16(), std::exception);
    EXPECT_THROW(stringTag.ToInt16(), std::exception);
    EXPECT_THROW(arrayTag.ToInt16(), std::exception);
}

TEST(TagTest, ConvertToUInt16)
{
    Tag zero(int32_t{0});
    Tag neg(int32_t{-1});
    Tag normal(int32_t{42});
    Tag max(int32_t{std::numeric_limits<uint16_t>::max()});
    Tag floatTag(float{3.14});
    Tag stringTag(std::string{"foo"});
    Tag arrayTag(std::vector<uint16_t>{{1, 2, 3}});

    // allowed
    EXPECT_NO_THROW({
        zero.ToUInt16();
        normal.ToUInt16();
        max.ToUInt16();
    });

    // not allowed
    EXPECT_THROW(neg.ToUInt16(), std::exception);
    EXPECT_THROW(floatTag.ToUInt16(), std::exception);
    EXPECT_THROW(stringTag.ToUInt16(), std::exception);
    EXPECT_THROW(arrayTag.ToUInt16(), std::exception);
}

TEST(TagTest, ConvertToInt32)
{
    Tag zero(int32_t{0});
    Tag min(int32_t{std::numeric_limits<int32_t>::min()});
    Tag normal(int32_t{42});
    Tag max(int32_t{std::numeric_limits<int32_t>::max()});
    Tag floatTag(float{3.14});
    Tag stringTag(std::string{"foo"});
    Tag arrayTag(std::vector<int32_t>{{1, 2, 3}});

    // allowed
    EXPECT_NO_THROW({
        zero.ToInt32();
        min.ToInt32();
        normal.ToInt32();
        max.ToInt32();
    });

    // not allowed
    EXPECT_THROW(floatTag.ToInt32(), std::exception);
    EXPECT_THROW(stringTag.ToInt32(), std::exception);
    EXPECT_THROW(arrayTag.ToInt32(), std::exception);
}

TEST(TagTest, ConvertToUInt32)
{
    Tag zero(int32_t{0});
    Tag neg(int32_t{-1});
    Tag normal(int32_t{42});
    Tag max(uint32_t{std::numeric_limits<uint32_t>::max()});
    Tag floatTag(float{3.14});
    Tag stringTag(std::string{"foo"});
    Tag arrayTag(std::vector<uint32_t>{{1, 2, 3}});

    // allowed
    EXPECT_NO_THROW({
        zero.ToUInt32();
        normal.ToUInt32();
        max.ToUInt32();
    });

    // not allowed
    EXPECT_THROW(neg.ToUInt32(), std::exception);
    EXPECT_THROW(floatTag.ToUInt32(), std::exception);
    EXPECT_THROW(stringTag.ToUInt32(), std::exception);
    EXPECT_THROW(arrayTag.ToUInt32(), std::exception);
}

TEST(TagCollectionTest, DefaultConstruction)
{
    TagCollection tags;
    EXPECT_TRUE(tags.empty());
    EXPECT_FALSE(tags.Contains("XY"));
}

TEST(TagCollectionTest, AddSimpleTags)
{
    const int32_t intValue = -42;
    const std::string strValue = "foo";
    const std::string hexStrValue = "1abc75";

    TagCollection tags;
    tags["ST"] = strValue;
    tags["XY"] = intValue;
    tags["HX"] = hexStrValue;
    tags["HX"].Modifier(TagModifier::HEX_STRING);

    EXPECT_EQ(3, tags.size());
    EXPECT_TRUE(tags.Contains("XY"));
    EXPECT_TRUE(tags.Contains("ST"));
    EXPECT_TRUE(tags.Contains("HX"));
    EXPECT_FALSE(tags.Contains("ZZ"));

    EXPECT_TRUE(tags["XY"].ToInt32() == intValue);
    EXPECT_TRUE(tags["ST"].ToString() == strValue);
    EXPECT_TRUE(tags["HX"].ToString() == hexStrValue);
    EXPECT_TRUE(tags["HX"].HasModifier(TagModifier::HEX_STRING));
}

TEST(SamTagCodecTest, DecodeTest)
{
    std::string tagString;
    tagString.append("HX:H:1abc75");
    tagString.append("\t");
    tagString.append("ST:Z:foo");
    tagString.append("\t");
    tagString.append("VC:B:i,42,-100,37,2048");
    tagString.append("\t");
    tagString.append("XY:i:-42");

    TagCollection expected;
    expected["ST"] = std::string("foo");
    expected["XY"] = int32_t{-42};
    expected["HX"] = Tag("1abc75", TagModifier::HEX_STRING);
    expected["VC"] = std::vector<int32_t>({42, -100, 37, 2048});

    TagCollection tags = SamTagCodec::Decode(tagString);

    EXPECT_TRUE(tags.Contains("ST"));
    EXPECT_TRUE(tags.Contains("HX"));
    EXPECT_TRUE(tags.Contains("XY"));
    EXPECT_TRUE(tags.Contains("VC"));

    EXPECT_EQ(std::string("foo"), tags["ST"].ToString());
    EXPECT_TRUE(tags["HX"].HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), tags["HX"].ToString());
    EXPECT_EQ(int8_t{-42}, tags["XY"].ToInt8());
    EXPECT_EQ(std::vector<int32_t>({42, -100, 37, 2048}), tags["VC"].ToInt32Array());
}

TEST(SamTagCodecTest, EncodeTest)
{
    TagCollection tags;
    tags["ST"] = std::string("foo");
    tags["XY"] = int32_t{-42};
    tags["HX"] = Tag("1abc75", TagModifier::HEX_STRING);
    tags["VC"] = std::vector<int32_t>({42, -100, 37, 2048});

    // "HX:H:1abc75\tST:Z:foo\0\tVC:B:i,42,-100,37,2048\tXY:i:-42"
    std::string expected;
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("ST:Z:foo");
    expected.append("\t");
    expected.append("VC:B:i,42,-100,37,2048");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(tags);
    EXPECT_EQ(expected, sam);
}

TEST(BamTagCodecTest, DecodeTest)
{
    std::vector<uint8_t> data;
    data.push_back(uint8_t('H'));
    data.push_back(uint8_t('X'));
    data.push_back(uint8_t('H'));
    data.push_back(uint8_t('1'));
    data.push_back(uint8_t('a'));
    data.push_back(uint8_t('b'));
    data.push_back(uint8_t('c'));
    data.push_back(uint8_t('7'));
    data.push_back(uint8_t('5'));
    data.push_back(uint8_t(0));

    data.push_back(uint8_t('X'));
    data.push_back(uint8_t('Y'));
    data.push_back(uint8_t('i'));
    const int32_t x = -42;
    char valueBytes[sizeof x];
    std::copy(static_cast<const char*>(static_cast<const void*>(&x)),
              static_cast<const char*>(static_cast<const void*>(&x)) + sizeof x, valueBytes);
    data.push_back(valueBytes[0]);
    data.push_back(valueBytes[1]);
    data.push_back(valueBytes[2]);
    data.push_back(valueBytes[3]);

    data.push_back('C');
    data.push_back('A');
    data.push_back('B');
    data.push_back('C');
    const uint32_t numChars = 3;
    char numCharsValueBytes[sizeof numChars];
    std::copy(static_cast<const char*>(static_cast<const void*>(&numChars)),
              static_cast<const char*>(static_cast<const void*>(&numChars)) + sizeof numChars,
              numCharsValueBytes);
    data.push_back(numCharsValueBytes[0]);
    data.push_back(numCharsValueBytes[1]);
    data.push_back(numCharsValueBytes[2]);
    data.push_back(numCharsValueBytes[3]);

    const std::vector<uint8_t> charArray = std::vector<uint8_t>({34, 5, 125});
    data.push_back(charArray.at(0));
    data.push_back(charArray.at(1));
    data.push_back(charArray.at(2));

    TagCollection tags = BamTagCodec::Decode(data);

    EXPECT_TRUE(tags["HX"].HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), tags["HX"].ToString());
    EXPECT_EQ(x, tags["XY"].ToInt32());
    EXPECT_EQ(charArray, tags["CA"].ToUInt8Array());

    // sanity check - convert tags back to SAM
    std::string expected;
    expected.append("CA:B:C,34,5,125");
    expected.append("\t");
    expected.append("HX:H:1abc75");
    expected.append("\t");
    expected.append("XY:i:-42");

    const std::string sam = SamTagCodec::Encode(tags);
    EXPECT_EQ(expected, sam);
}

TEST(BamTagCodecTest, EncodeTest)
{
    std::vector<uint8_t> expected;

    expected.push_back('C');
    expected.push_back('A');
    expected.push_back('B');
    expected.push_back('C');
    const uint32_t numChars = 3;
    char numCharsValueBytes[sizeof numChars];
    std::copy(static_cast<const char*>(static_cast<const void*>(&numChars)),
              static_cast<const char*>(static_cast<const void*>(&numChars)) + sizeof numChars,
              numCharsValueBytes);
    expected.push_back(numCharsValueBytes[0]);
    expected.push_back(numCharsValueBytes[1]);
    expected.push_back(numCharsValueBytes[2]);
    expected.push_back(numCharsValueBytes[3]);

    const std::vector<uint8_t> charArray = std::vector<uint8_t>({34, 5, 125});
    expected.push_back(charArray.at(0));
    expected.push_back(charArray.at(1));
    expected.push_back(charArray.at(2));

    expected.push_back(uint8_t('H'));
    expected.push_back(uint8_t('X'));
    expected.push_back(uint8_t('H'));
    expected.push_back(uint8_t('1'));
    expected.push_back(uint8_t('a'));
    expected.push_back(uint8_t('b'));
    expected.push_back(uint8_t('c'));
    expected.push_back(uint8_t('7'));
    expected.push_back(uint8_t('5'));
    expected.push_back(uint8_t(0));

    expected.push_back(uint8_t('X'));
    expected.push_back(uint8_t('Y'));
    expected.push_back(uint8_t('i'));
    const int32_t x = -42;
    char valueBytes[sizeof x];
    std::copy(static_cast<const char*>(static_cast<const void*>(&x)),
              static_cast<const char*>(static_cast<const void*>(&x)) + sizeof x, valueBytes);
    expected.push_back(valueBytes[0]);
    expected.push_back(valueBytes[1]);
    expected.push_back(valueBytes[2]);
    expected.push_back(valueBytes[3]);

    TagCollection tags;
    tags["HX"] = Tag("1abc75", TagModifier::HEX_STRING);
    tags["CA"] = charArray;
    tags["XY"] = x;

    const std::vector<uint8_t> data = BamTagCodec::Encode(tags);
    EXPECT_EQ(expected, data);
}

TEST(BamTagCodecTest, AsciiTagsTest)
{
    std::vector<uint8_t> expected;
    expected.reserve(20);
    expected.push_back('I');  // I8:A:B
    expected.push_back('8');
    expected.push_back('A');
    expected.push_back('B');
    expected.push_back('P');  // PC:A:$
    expected.push_back('C');
    expected.push_back('A');
    expected.push_back('$');
    expected.push_back('S');  // SC:A:$
    expected.push_back('C');
    expected.push_back('A');
    expected.push_back('$');
    expected.push_back('U');  // U8:A:A
    expected.push_back('8');
    expected.push_back('A');
    expected.push_back('A');
    expected.push_back('U');  // UC:A:$
    expected.push_back('C');
    expected.push_back('A');
    expected.push_back('$');

    const char c = '$';
    const signed char sc = '$';
    const unsigned char uc = '$';
    const uint8_t u8 = 65;
    const int8_t i8 = 66;

    {  // old style: construct-then-modify

        Tag fromPlainChar = Tag(c);
        Tag fromSignedChar = Tag(sc);
        Tag fromUnsignedChar = Tag(uc);
        Tag fromUint8 = Tag(u8);
        Tag fromInt8 = Tag(i8);
        fromPlainChar.Modifier(TagModifier::ASCII_CHAR);
        fromSignedChar.Modifier(TagModifier::ASCII_CHAR);
        fromUnsignedChar.Modifier(TagModifier::ASCII_CHAR);
        fromUint8.Modifier(TagModifier::ASCII_CHAR);
        fromInt8.Modifier(TagModifier::ASCII_CHAR);

        TagCollection tags;
        tags["PC"] = fromPlainChar;
        tags["SC"] = fromSignedChar;
        tags["UC"] = fromUnsignedChar;
        tags["U8"] = fromUint8;
        tags["I8"] = fromInt8;

        const std::vector<uint8_t> data = BamTagCodec::Encode(tags);
        EXPECT_EQ(expected, data);
    }

    {  // new style: construct directly as ASCII

        const Tag fromPlainChar = Tag(c, TagModifier::ASCII_CHAR);
        const Tag fromSignedChar = Tag(sc, TagModifier::ASCII_CHAR);
        const Tag fromUnsignedChar = Tag(uc, TagModifier::ASCII_CHAR);
        const Tag fromUint8 = Tag(u8, TagModifier::ASCII_CHAR);
        const Tag fromInt8 = Tag(i8, TagModifier::ASCII_CHAR);

        TagCollection tags;
        tags["PC"] = fromPlainChar;
        tags["SC"] = fromSignedChar;
        tags["UC"] = fromUnsignedChar;
        tags["U8"] = fromUint8;
        tags["I8"] = fromInt8;

        const std::vector<uint8_t> data = BamTagCodec::Encode(tags);
        EXPECT_EQ(expected, data);
    }
}

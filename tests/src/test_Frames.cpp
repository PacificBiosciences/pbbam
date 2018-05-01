// Author: Derek Barnett

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/Frames.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace FramesTests {

static const std::vector<uint16_t> testFrames{
    0,  8,  140, 0,  0,   7,  4,  0,  85, 2,  1,  3,  2,   10, 1,  20, 47,   10,  9,  60, 20,
    3,  12, 5,   13, 165, 6,  14, 22, 12, 2,  4,  9,  218, 27, 3,  15, 2,    17,  2,  45, 24,
    89, 10, 7,   1,  11,  15, 0,  7,  0,  28, 17, 12, 6,   10, 37, 0,  12,   52,  0,  7,  1,
    14, 3,  26,  12, 0,   20, 17, 2,  13, 2,  9,  13, 7,   15, 29, 3,  6,    2,   1,  28, 10,
    3,  14, 7,   1,  22,  1,  6,  6,  0,  19, 31, 6,  2,   14, 0,  0,  1000, 947, 948};

static const std::vector<uint8_t> encodedFrames{
    0,  8,  102, 0,  0,   7,  4,  0,  75, 2,  1,  3,  2,   10, 1,  20, 47,  10,  9,  60, 20,
    3,  12, 5,   13, 115, 6,  14, 22, 12, 2,  4,  9,  135, 27, 3,  15, 2,   17,  2,  45, 24,
    77, 10, 7,   1,  11,  15, 0,  7,  0,  28, 17, 12, 6,   10, 37, 0,  12,  52,  0,  7,  1,
    14, 3,  26,  12, 0,   20, 17, 2,  13, 2,  9,  13, 7,   15, 29, 3,  6,   2,   1,  28, 10,
    3,  14, 7,   1,  22,  1,  6,  6,  0,  19, 31, 6,  2,   14, 0,  0,  255, 254, 255};

}  // namespace FramesTests

TEST(FramesTest, Constructors)
{
    const Frames f;
    ASSERT_TRUE(f.Data().empty());

    const Frames f2(FramesTests::testFrames);
    const auto d = f2.Data();
    ASSERT_EQ(FramesTests::testFrames, d);
}

TEST(FramesTest, Encoded)
{
    const Frames f(FramesTests::testFrames);
    const auto e = f.Encode();
    ASSERT_EQ(FramesTests::encodedFrames, e);
}

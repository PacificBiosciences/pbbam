// Author: Derek Barnett

#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/Compare.h>

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;

namespace CompareTests {

static inline
BamRecord makeRecordWithTag(const std::string& tagName,
                            const Tag& tag)
{
    auto r = BamRecord{ };
    r.Impl().AddTag(tagName, tag);
    return r;
}

static
BamRecord makeRecord(const Position qStart,
                     const Position qEnd,
                     const std::string& seq,
                     const std::string& quals,
                     const std::string& tagBases,
                     const std::string& tagQuals,
                     const std::vector<uint16_t>& frames)
{
    BamRecordImpl impl;
    impl.SetSequenceAndQualities(seq, quals);

    TagCollection tags;
    tags["qs"] = qStart;
    tags["qe"] = qEnd;
    tags["ip"] = frames;
    tags["pw"] = frames;
    tags["dt"] = tagBases;
    tags["st"] = tagBases;
    tags["dq"] = tagQuals;
    tags["iq"] = tagQuals;
    tags["mq"] = tagQuals;
    tags["sq"] = tagQuals;
    tags["pq"] = tagQuals;
    tags["pv"] = tagQuals;
    impl.Tags(tags);

    return BamRecord(std::move(impl));
}

static
std::vector<BamRecord> makeMappedRecords()
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const std::vector<uint16_t> frames  = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const uint8_t mapQual = 80;

    const std::string s1_cigar = "10=";
    const std::string s2_cigar = "5=3D5=";
    const std::string s3_cigar = "4=1D2I2D2X2=";

    BamRecord s1 = CompareTests::makeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s2 = CompareTests::makeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s3 = CompareTests::makeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s1_rev = CompareTests::makeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s2_rev = CompareTests::makeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s3_rev = CompareTests::makeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);

    s1.Map(0, 100, Strand::FORWARD, s1_cigar, mapQual);
    s2.Map(0, 100, Strand::FORWARD, s2_cigar, mapQual);
    s3.Map(0, 100, Strand::FORWARD, s3_cigar, mapQual);
    s1_rev.Map(0, 100, Strand::REVERSE, s1_cigar, mapQual);
    s2_rev.Map(0, 100, Strand::REVERSE, s2_cigar, mapQual);
    s3_rev.Map(0, 100, Strand::REVERSE, s3_cigar, mapQual);

    return std::vector<BamRecord> { s1, s2, s3, s1_rev, s2_rev, s3_rev };
}

} // namespace CompareTests

TEST(CompareTest, TypeToNameOk)
{
    EXPECT_EQ(std::string{"Compare::EQUAL"},              Compare::TypeToName(Compare::EQUAL));
    EXPECT_EQ(std::string{"Compare::NOT_EQUAL"},          Compare::TypeToName(Compare::NOT_EQUAL));
    EXPECT_EQ(std::string{"Compare::LESS_THAN"},          Compare::TypeToName(Compare::LESS_THAN));
    EXPECT_EQ(std::string{"Compare::LESS_THAN_EQUAL"},    Compare::TypeToName(Compare::LESS_THAN_EQUAL));
    EXPECT_EQ(std::string{"Compare::GREATER_THAN"},       Compare::TypeToName(Compare::GREATER_THAN));
    EXPECT_EQ(std::string{"Compare::GREATER_THAN_EQUAL"}, Compare::TypeToName(Compare::GREATER_THAN_EQUAL));
    EXPECT_EQ(std::string{"Compare::CONTAINS"},           Compare::TypeToName(Compare::CONTAINS));
    EXPECT_EQ(std::string{"Compare::NOT_CONTAINS"},       Compare::TypeToName(Compare::NOT_CONTAINS));

    // invalid type throws
    EXPECT_THROW(Compare::TypeToName(static_cast<Compare::Type>(42)), std::runtime_error);
}

TEST(CompareTest, TypeToOperatorOk)
{
    { // normal
        EXPECT_EQ(Compare::TypeToOperator(Compare::EQUAL),              std::string{"=="});
        EXPECT_EQ(Compare::TypeToOperator(Compare::NOT_EQUAL),          std::string{"!="});
        EXPECT_EQ(Compare::TypeToOperator(Compare::LESS_THAN),          std::string{"<"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::LESS_THAN_EQUAL),    std::string{"<="});
        EXPECT_EQ(Compare::TypeToOperator(Compare::GREATER_THAN),       std::string{">"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::GREATER_THAN_EQUAL), std::string{">="});
        EXPECT_EQ(Compare::TypeToOperator(Compare::CONTAINS),           std::string{"&"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::NOT_CONTAINS),       std::string{"~"});
    }

    { // alpha
        EXPECT_EQ(Compare::TypeToOperator(Compare::EQUAL, true),              std::string{"eq"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::NOT_EQUAL, true),          std::string{"ne"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::LESS_THAN, true),          std::string{"lt"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::LESS_THAN_EQUAL, true),    std::string{"lte"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::GREATER_THAN, true),       std::string{"gt"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::GREATER_THAN_EQUAL, true), std::string{"gte"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::CONTAINS, true),           std::string{"and"});
        EXPECT_EQ(Compare::TypeToOperator(Compare::NOT_CONTAINS, true),       std::string{"not"});
    }

    // invalid type throws
    EXPECT_THROW(Compare::TypeToOperator(static_cast<Compare::Type>(42)), std::runtime_error);
}

TEST(CompareTest, FromOperatorOk)
{
    EXPECT_EQ(Compare::EQUAL,              Compare::TypeFromOperator("=="));
    EXPECT_EQ(Compare::EQUAL,              Compare::TypeFromOperator("="));
    EXPECT_EQ(Compare::EQUAL,              Compare::TypeFromOperator("eq"));
    EXPECT_EQ(Compare::NOT_EQUAL,          Compare::TypeFromOperator("!="));
    EXPECT_EQ(Compare::NOT_EQUAL,          Compare::TypeFromOperator("ne"));
    EXPECT_EQ(Compare::LESS_THAN,          Compare::TypeFromOperator("<"));
    EXPECT_EQ(Compare::LESS_THAN,          Compare::TypeFromOperator("lt"));
    EXPECT_EQ(Compare::LESS_THAN,          Compare::TypeFromOperator("&lt;"));
    EXPECT_EQ(Compare::LESS_THAN_EQUAL,    Compare::TypeFromOperator("<="));
    EXPECT_EQ(Compare::LESS_THAN_EQUAL,    Compare::TypeFromOperator("lte"));
    EXPECT_EQ(Compare::LESS_THAN_EQUAL,    Compare::TypeFromOperator("&lt;="));
    EXPECT_EQ(Compare::GREATER_THAN,       Compare::TypeFromOperator(">"));
    EXPECT_EQ(Compare::GREATER_THAN,       Compare::TypeFromOperator("gt"));
    EXPECT_EQ(Compare::GREATER_THAN,       Compare::TypeFromOperator("&gt;"));
    EXPECT_EQ(Compare::GREATER_THAN_EQUAL, Compare::TypeFromOperator(">="));
    EXPECT_EQ(Compare::GREATER_THAN_EQUAL, Compare::TypeFromOperator("gte"));
    EXPECT_EQ(Compare::GREATER_THAN_EQUAL, Compare::TypeFromOperator("&gt;="));
    EXPECT_EQ(Compare::CONTAINS,           Compare::TypeFromOperator("&"));
    EXPECT_EQ(Compare::NOT_CONTAINS,       Compare::TypeFromOperator("~"));

    // invalid operator strings throw
    EXPECT_THROW(Compare::TypeFromOperator(""),        std::runtime_error);
    EXPECT_THROW(Compare::TypeFromOperator("invalid"), std::runtime_error);
}

TEST(CompareTest, AlignedEndOk)
{
    BamRecord r1;
    r1.Map(0, 290, Strand::FORWARD, Cigar{"10="}, 255);

    BamRecord r2;
    r2.Map(0, 190, Strand::FORWARD, Cigar{"10="}, 255);

    BamRecord r3;
    r3.Map(0, 290, Strand::FORWARD, Cigar{"10="}, 255);

    BamRecord r4;
    r4.Map(0, 90, Strand::FORWARD, Cigar{"10="}, 255);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::AlignedEnd());

    EXPECT_EQ(r4.AlignedEnd(), records.at(0).AlignedEnd());
    EXPECT_EQ(r2.AlignedEnd(), records.at(1).AlignedEnd());
    EXPECT_EQ(r1.AlignedEnd(), records.at(2).AlignedEnd());
    EXPECT_EQ(r3.AlignedEnd(), records.at(3).AlignedEnd());
}

TEST(CompareTest, AlignedStartOk)
{
    BamRecord r1;
    r1.Map(0, 300, Strand::FORWARD, Cigar{"10="}, 255);

    BamRecord r2;
    r2.Map(0, 200, Strand::FORWARD, Cigar{"10="}, 255);

    BamRecord r3;
    r3.Map(0, 400, Strand::FORWARD, Cigar{"10="}, 255);

    BamRecord r4;
    r4.Map(0, 100, Strand::FORWARD, Cigar{"10="}, 255);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::AlignedStart());

    EXPECT_EQ(r4.AlignedStart(), records.at(0).AlignedStart());
    EXPECT_EQ(r2.AlignedStart(), records.at(1).AlignedStart());
    EXPECT_EQ(r1.AlignedStart(), records.at(2).AlignedStart());
    EXPECT_EQ(r3.AlignedStart(), records.at(3).AlignedStart());
}

TEST(CompareTest, AlignedStrandOk)
{
    BamRecord r1; r1.Impl().SetReverseStrand(true);
    BamRecord r2; r2.Impl().SetReverseStrand(false);
    BamRecord r3; r3.Impl().SetReverseStrand(true);
    BamRecord r4; r4.Impl().SetReverseStrand(false);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::AlignedStrand());

    EXPECT_EQ(Strand::FORWARD, records.at(0).AlignedStrand());
    EXPECT_EQ(Strand::FORWARD, records.at(1).AlignedStrand());
    EXPECT_EQ(Strand::REVERSE, records.at(2).AlignedStrand());
    EXPECT_EQ(Strand::REVERSE, records.at(3).AlignedStrand());
}

TEST(CompareTest, BarcodeForwardOk)
{
    BamRecord r1; r1.Barcodes(std::make_pair<int16_t,int16_t>(30,20));
    BamRecord r2; r2.Barcodes(std::make_pair<int16_t,int16_t>(20,30));
    BamRecord r3; r3.Barcodes(std::make_pair<int16_t,int16_t>(40,10));
    BamRecord r4; r4.Barcodes(std::make_pair<int16_t,int16_t>(10,40));

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::BarcodeForward());

    EXPECT_EQ(r4.BarcodeForward(), records.at(0).BarcodeForward());
    EXPECT_EQ(r2.BarcodeForward(), records.at(1).BarcodeForward());
    EXPECT_EQ(r1.BarcodeForward(), records.at(2).BarcodeForward());
    EXPECT_EQ(r3.BarcodeForward(), records.at(3).BarcodeForward());
}

TEST(CompareTest, BarcodeReverseOk)
{
    BamRecord r1; r1.Barcodes(std::make_pair<int16_t,int16_t>(30,20));
    BamRecord r2; r2.Barcodes(std::make_pair<int16_t,int16_t>(20,30));
    BamRecord r3; r3.Barcodes(std::make_pair<int16_t,int16_t>(40,10));
    BamRecord r4; r4.Barcodes(std::make_pair<int16_t,int16_t>(10,40));

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::BarcodeReverse());

    EXPECT_EQ(r3.BarcodeReverse(), records.at(0).BarcodeReverse());
    EXPECT_EQ(r1.BarcodeReverse(), records.at(1).BarcodeReverse());
    EXPECT_EQ(r2.BarcodeReverse(), records.at(2).BarcodeReverse());
    EXPECT_EQ(r4.BarcodeReverse(), records.at(3).BarcodeReverse());
}

TEST(CompareTest, BarcodeQualityOk)
{
    uint8_t q1 = 30;
    uint8_t q2 = 20;
    uint8_t q3 = 40;
    uint8_t q4 = 10;

    auto records = std::vector<BamRecord>
    {
        CompareTests::makeRecordWithTag("bq", Tag(q1)),
        CompareTests::makeRecordWithTag("bq", Tag(q2)),
        CompareTests::makeRecordWithTag("bq", Tag(q3)),
        CompareTests::makeRecordWithTag("bq", Tag(q4))
    };
    std::sort(records.begin(), records.end(), Compare::BarcodeQuality());

    EXPECT_EQ(q4, records.at(0).BarcodeQuality());
    EXPECT_EQ(q2, records.at(1).BarcodeQuality());
    EXPECT_EQ(q1, records.at(2).BarcodeQuality());
    EXPECT_EQ(q3, records.at(3).BarcodeQuality());
}

TEST(CompareTest, CustomCompareOk)
{
    struct CustomCompare : public Compare::MemberFunctionBase<bool, &BamRecord::HasDeletionTag> { };

    auto records = std::vector<BamRecord>
    {
        CompareTests::makeRecordWithTag("dt", Tag(std::string("foo"))),
        CompareTests::makeRecordWithTag("dt", Tag(std::string("foo"))),
        CompareTests::makeRecordWithTag("dt", Tag(std::string("foo"))),
        CompareTests::makeRecordWithTag("dt", Tag(std::string("foo")))
    };
    records.push_back(BamRecord());
    records.push_back(BamRecord());
    records.push_back(BamRecord());
    records.push_back(BamRecord());
    EXPECT_EQ(8, records.size());

    std::sort(records.begin(), records.end(), CustomCompare());

    EXPECT_FALSE(records.at(0).HasDeletionTag());
    EXPECT_FALSE(records.at(1).HasDeletionTag());
    EXPECT_FALSE(records.at(2).HasDeletionTag());
    EXPECT_FALSE(records.at(3).HasDeletionTag());
    EXPECT_TRUE(records.at(4).HasDeletionTag());
    EXPECT_TRUE(records.at(5).HasDeletionTag());
    EXPECT_TRUE(records.at(6).HasDeletionTag());
    EXPECT_TRUE(records.at(7).HasDeletionTag());
}

TEST(CompareTest, FullNameOk)
{
    BamRecord r1; r1.Impl().Name("c");
    BamRecord r2; r2.Impl().Name("b");
    BamRecord r3; r3.Impl().Name("d");
    BamRecord r4; r4.Impl().Name("a");

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::FullName());

    EXPECT_EQ(r4.FullName(), records.at(0).FullName());
    EXPECT_EQ(r2.FullName(), records.at(1).FullName());
    EXPECT_EQ(r1.FullName(), records.at(2).FullName());
    EXPECT_EQ(r3.FullName(), records.at(3).FullName());
}

TEST(CompareTest, LocalContextFlagOk)
{
    BamRecord r1; r1.LocalContextFlags(LocalContextFlags::BARCODE_AFTER);
    BamRecord r2; r2.LocalContextFlags(LocalContextFlags::ADAPTER_AFTER);
    BamRecord r3; r3.LocalContextFlags(LocalContextFlags::REVERSE_PASS);
    BamRecord r4; r4.LocalContextFlags(LocalContextFlags::NO_LOCAL_CONTEXT);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::LocalContextFlag());

    EXPECT_EQ(r4.LocalContextFlags(), records.at(0).LocalContextFlags());
    EXPECT_EQ(r2.LocalContextFlags(), records.at(1).LocalContextFlags());
    EXPECT_EQ(r1.LocalContextFlags(), records.at(2).LocalContextFlags());
    EXPECT_EQ(r3.LocalContextFlags(), records.at(3).LocalContextFlags());
}

TEST(CompareTest, MapQualityOk)
{
    BamRecord r1; r1.Impl().MapQuality(30);
    BamRecord r2; r2.Impl().MapQuality(20);
    BamRecord r3; r3.Impl().MapQuality(40);
    BamRecord r4; r4.Impl().MapQuality(10);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::MapQuality());

    EXPECT_EQ(r4.MapQuality(), records.at(0).MapQuality());
    EXPECT_EQ(r2.MapQuality(), records.at(1).MapQuality());
    EXPECT_EQ(r1.MapQuality(), records.at(2).MapQuality());
    EXPECT_EQ(r3.MapQuality(), records.at(3).MapQuality());
}

TEST(CompareTest, MovieNameOk)
{
    auto rg1 = ReadGroupInfo { "a", "SUBREAD" };
    auto rg2 = ReadGroupInfo { "b", "SUBREAD" };
    auto rg3 = ReadGroupInfo { "c", "SUBREAD" };
    auto rg4 = ReadGroupInfo { "d", "SUBREAD" };

    BamHeader header;
    header.AddReadGroup(rg1)
          .AddReadGroup(rg2)
          .AddReadGroup(rg3)
          .AddReadGroup(rg4);

    BamRecord r1(header); r1.ReadGroup(rg3);
    BamRecord r2(header); r2.ReadGroup(rg2);
    BamRecord r3(header); r3.ReadGroup(rg4);
    BamRecord r4(header); r4.ReadGroup(rg1);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::MovieName());

    EXPECT_EQ(r4.MovieName(), records.at(0).MovieName());
    EXPECT_EQ(r2.MovieName(), records.at(1).MovieName());
    EXPECT_EQ(r1.MovieName(), records.at(2).MovieName());
    EXPECT_EQ(r3.MovieName(), records.at(3).MovieName());
}

TEST(CompareTest, NoneOk)
{
    BamRecord r1; r1.Impl().Name("c");
    BamRecord r2; r2.Impl().Name("b");
    BamRecord r3; r3.Impl().Name("d");
    BamRecord r4; r4.Impl().Name("a");

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::None());

    EXPECT_EQ(r1.FullName(), records.at(0).FullName());
    EXPECT_EQ(r2.FullName(), records.at(1).FullName());
    EXPECT_EQ(r3.FullName(), records.at(2).FullName());
    EXPECT_EQ(r4.FullName(), records.at(3).FullName());
}

TEST(CompareTest, NumDeletedBasesOk)
{
    // create test data
    auto records = CompareTests::makeMappedRecords();

    // sanity checks on initial conditions
    EXPECT_EQ(6, records.size());
    EXPECT_EQ(0, records.at(0).NumDeletedBases());
    EXPECT_EQ(3, records.at(1).NumDeletedBases());
    EXPECT_EQ(3, records.at(2).NumDeletedBases());
    EXPECT_EQ(0, records.at(3).NumDeletedBases());
    EXPECT_EQ(3, records.at(4).NumDeletedBases());
    EXPECT_EQ(3, records.at(5).NumDeletedBases());

    // sort & check
    std::sort(records.begin(), records.end(), Compare::NumDeletedBases());
    EXPECT_EQ(0, records.at(0).NumDeletedBases());
    EXPECT_EQ(0, records.at(1).NumDeletedBases());
    EXPECT_EQ(3, records.at(2).NumDeletedBases());
    EXPECT_EQ(3, records.at(3).NumDeletedBases());
    EXPECT_EQ(3, records.at(4).NumDeletedBases());
    EXPECT_EQ(3, records.at(5).NumDeletedBases());
}

TEST(CompareTest, NumInsertedBasesOk)
{
    // create test data
    auto records = CompareTests::makeMappedRecords();

    // sanity checks on initial conditions
    EXPECT_EQ(6, records.size());
    EXPECT_EQ(0, records.at(0).NumInsertedBases());
    EXPECT_EQ(0, records.at(1).NumInsertedBases());
    EXPECT_EQ(2, records.at(2).NumInsertedBases());
    EXPECT_EQ(0, records.at(3).NumInsertedBases());
    EXPECT_EQ(0, records.at(4).NumInsertedBases());
    EXPECT_EQ(2, records.at(5).NumInsertedBases());

    // sort & check
    std::sort(records.begin(), records.end(), Compare::NumInsertedBases());
    EXPECT_EQ(0, records.at(0).NumInsertedBases());
    EXPECT_EQ(0, records.at(1).NumInsertedBases());
    EXPECT_EQ(0, records.at(2).NumInsertedBases());
    EXPECT_EQ(0, records.at(3).NumInsertedBases());
    EXPECT_EQ(2, records.at(4).NumInsertedBases());
    EXPECT_EQ(2, records.at(5).NumInsertedBases());
}

TEST(CompareTest, NumMatchesOk)
{
    // create test data
    auto records = CompareTests::makeMappedRecords();

    // sanity checks on initial conditions
    EXPECT_EQ(6, records.size());
    EXPECT_EQ(10, records.at(0).NumMatches());
    EXPECT_EQ(10, records.at(1).NumMatches());
    EXPECT_EQ(6,  records.at(2).NumMatches());
    EXPECT_EQ(10, records.at(3).NumMatches());
    EXPECT_EQ(10, records.at(4).NumMatches());
    EXPECT_EQ(6,  records.at(5).NumMatches());

    // sort & check
    std::sort(records.begin(), records.end(), Compare::NumMatches());
    EXPECT_EQ(6,  records.at(0).NumMatches());
    EXPECT_EQ(6,  records.at(1).NumMatches());
    EXPECT_EQ(10, records.at(2).NumMatches());
    EXPECT_EQ(10, records.at(3).NumMatches());
    EXPECT_EQ(10, records.at(4).NumMatches());
    EXPECT_EQ(10, records.at(5).NumMatches());
}

TEST(CompareTest, NumMismatchesOk)
{
    // create test data
    auto records = CompareTests::makeMappedRecords();

    // sanity checks on initial conditions
    EXPECT_EQ(6, records.size());
    EXPECT_EQ(0, records.at(0).NumMismatches());
    EXPECT_EQ(0, records.at(1).NumMismatches());
    EXPECT_EQ(2, records.at(2).NumMismatches());
    EXPECT_EQ(0, records.at(3).NumMismatches());
    EXPECT_EQ(0, records.at(4).NumMismatches());
    EXPECT_EQ(2, records.at(5).NumMismatches());

    // sort & check
    std::sort(records.begin(), records.end(), Compare::NumMismatches());
    EXPECT_EQ(0, records.at(0).NumMismatches());
    EXPECT_EQ(0, records.at(1).NumMismatches());
    EXPECT_EQ(0, records.at(2).NumMismatches());
    EXPECT_EQ(0, records.at(3).NumMismatches());
    EXPECT_EQ(2, records.at(4).NumMismatches());
    EXPECT_EQ(2, records.at(5).NumMismatches());
}

TEST(CompareTest, QueryEndOk)
{
    Position q1 = 30;
    Position q2 = 20;
    Position q3 = 40;
    Position q4 = 10;

    auto records = std::vector<BamRecord>
    {
        CompareTests::makeRecordWithTag("qe", Tag(q1)),
        CompareTests::makeRecordWithTag("qe", Tag(q2)),
        CompareTests::makeRecordWithTag("qe", Tag(q3)),
        CompareTests::makeRecordWithTag("qe", Tag(q4))
    };
    std::sort(records.begin(), records.end(), Compare::QueryEnd());

    EXPECT_EQ(q4, records.at(0).QueryEnd());
    EXPECT_EQ(q2, records.at(1).QueryEnd());
    EXPECT_EQ(q1, records.at(2).QueryEnd());
    EXPECT_EQ(q3, records.at(3).QueryEnd());
}

TEST(CompareTest, QueryStartOk)
{
    Position q1 = 30;
    Position q2 = 20;
    Position q3 = 40;
    Position q4 = 10;

    auto records = std::vector<BamRecord>
    {
        CompareTests::makeRecordWithTag("qs", Tag(q1)),
        CompareTests::makeRecordWithTag("qs", Tag(q2)),
        CompareTests::makeRecordWithTag("qs", Tag(q3)),
        CompareTests::makeRecordWithTag("qs", Tag(q4))
    };
    std::sort(records.begin(), records.end(), Compare::QueryStart());

    EXPECT_EQ(q4, records.at(0).QueryStart());
    EXPECT_EQ(q2, records.at(1).QueryStart());
    EXPECT_EQ(q1, records.at(2).QueryStart());
    EXPECT_EQ(q3, records.at(3).QueryStart());
}

TEST(CompareTest, ReadGroupIdOk)
{
    auto rg1 = ReadGroupInfo { "foo", "SUBREAD" };
    auto rg2 = ReadGroupInfo { "bar", "SUBREAD" };
    auto rg3 = ReadGroupInfo { "c",   "SUBREAD" };
    auto rg4 = ReadGroupInfo { "d",   "SUBREAD" };

    BamHeader header;
    header.AddReadGroup(rg1)
          .AddReadGroup(rg2)
          .AddReadGroup(rg3)
          .AddReadGroup(rg4);

    BamRecord r1(header); r1.ReadGroup(rg3); // -> 99365356
    BamRecord r2(header); r2.ReadGroup(rg2); // -> d9f305e4
    BamRecord r3(header); r3.ReadGroup(rg4); // -> 54397cd6
    BamRecord r4(header); r4.ReadGroup(rg1); // -> a60ddc69

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::ReadGroupId()); // lexical, NOT numeric ordering

    EXPECT_EQ(r3.ReadGroupId(), records.at(0).ReadGroupId());
    EXPECT_EQ(r1.ReadGroupId(), records.at(1).ReadGroupId());
    EXPECT_EQ(r4.ReadGroupId(), records.at(2).ReadGroupId());
    EXPECT_EQ(r2.ReadGroupId(), records.at(3).ReadGroupId());
}

TEST(CompareTest, ReadGroupNumericIdOk)
{
    auto rg1 = ReadGroupInfo { "a", "SUBREAD" };
    auto rg2 = ReadGroupInfo { "b", "SUBREAD" };
    auto rg3 = ReadGroupInfo { "c", "SUBREAD" };
    auto rg4 = ReadGroupInfo { "d", "SUBREAD" };

    BamHeader header;
    header.AddReadGroup(rg1)
          .AddReadGroup(rg2)
          .AddReadGroup(rg3)
          .AddReadGroup(rg4);

    BamRecord r1(header); r1.ReadGroup(rg3); // -> -1724492970
    BamRecord r2(header); r2.ReadGroup(rg2); // ->   235381373
    BamRecord r3(header); r3.ReadGroup(rg4); // ->  1413053654
    BamRecord r4(header); r4.ReadGroup(rg1); // ->  1153643386

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::ReadGroupNumericId()); // numeric ordering

    EXPECT_EQ(r1.ReadGroupNumericId(), records.at(0).ReadGroupNumericId());
    EXPECT_EQ(r2.ReadGroupNumericId(), records.at(1).ReadGroupNumericId());
    EXPECT_EQ(r4.ReadGroupNumericId(), records.at(2).ReadGroupNumericId());
    EXPECT_EQ(r3.ReadGroupNumericId(), records.at(3).ReadGroupNumericId());
}

TEST(CompareTest, ReadAccuracyOk)
{
    Accuracy a1 = 30;
    Accuracy a2 = 20;
    Accuracy a3 = 40;
    Accuracy a4 = 10;

    auto records = std::vector<BamRecord>
    {
        CompareTests::makeRecordWithTag("rq", Tag(a1)),
        CompareTests::makeRecordWithTag("rq", Tag(a2)),
        CompareTests::makeRecordWithTag("rq", Tag(a3)),
        CompareTests::makeRecordWithTag("rq", Tag(a4))
    };
    std::sort(records.begin(), records.end(), Compare::ReadAccuracy());

    EXPECT_EQ(a4, records.at(0).ReadAccuracy());
    EXPECT_EQ(a2, records.at(1).ReadAccuracy());
    EXPECT_EQ(a1, records.at(2).ReadAccuracy());
    EXPECT_EQ(a3, records.at(3).ReadAccuracy());
}

TEST(CompareTest, ReferenceEndOk)
{
    // create test data
    auto records = CompareTests::makeMappedRecords();

    // sanity checks on initial conditions
    EXPECT_EQ(6, records.size());
    EXPECT_EQ(110, records.at(0).ReferenceEnd());
    EXPECT_EQ(113, records.at(1).ReferenceEnd());
    EXPECT_EQ(111, records.at(2).ReferenceEnd());
    EXPECT_EQ(110, records.at(3).ReferenceEnd());
    EXPECT_EQ(113, records.at(4).ReferenceEnd());
    EXPECT_EQ(111, records.at(5).ReferenceEnd());

    // sort & check
    std::sort(records.begin(), records.end(), Compare::ReferenceEnd());
    EXPECT_EQ(110, records.at(0).ReferenceEnd());
    EXPECT_EQ(110, records.at(1).ReferenceEnd());
    EXPECT_EQ(111, records.at(2).ReferenceEnd());
    EXPECT_EQ(111, records.at(3).ReferenceEnd());
    EXPECT_EQ(113, records.at(4).ReferenceEnd());
    EXPECT_EQ(113, records.at(5).ReferenceEnd());
}

TEST(CompareTest, ReferenceIdOk)
{
    BamRecord r1; r1.Impl().ReferenceId(30);
    BamRecord r2; r2.Impl().ReferenceId(20);
    BamRecord r3; r3.Impl().ReferenceId(40);
    BamRecord r4; r4.Impl().ReferenceId(10);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::ReferenceId());

    EXPECT_EQ(r4.ReferenceId(), records.at(0).ReferenceId());
    EXPECT_EQ(r2.ReferenceId(), records.at(1).ReferenceId());
    EXPECT_EQ(r1.ReferenceId(), records.at(2).ReferenceId());
    EXPECT_EQ(r3.ReferenceId(), records.at(3).ReferenceId());
}

TEST(CompareTest, ReferenceNameOk)
{
    auto seq1 = SequenceInfo { "seq1" };
    auto seq2 = SequenceInfo { "seq2" };
    auto seq3 = SequenceInfo { "seq3" };
    auto seq4 = SequenceInfo { "seq4" };

    BamHeader header;
    header.AddSequence(seq1)  // -> 0
          .AddSequence(seq2)  // -> 1
          .AddSequence(seq3)  // -> 2
          .AddSequence(seq4); // -> 3

    BamRecord r1(header); r1.Impl().SetMapped(true); r1.Impl().ReferenceId(2);
    BamRecord r2(header); r2.Impl().SetMapped(true); r2.Impl().ReferenceId(1);
    BamRecord r3(header); r3.Impl().SetMapped(true); r3.Impl().ReferenceId(3);
    BamRecord r4(header); r4.Impl().SetMapped(true); r4.Impl().ReferenceId(0);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::ReferenceName());

    EXPECT_EQ(seq1.Name(), records.at(0).ReferenceName());
    EXPECT_EQ(seq2.Name(), records.at(1).ReferenceName());
    EXPECT_EQ(seq3.Name(), records.at(2).ReferenceName());
    EXPECT_EQ(seq4.Name(), records.at(3).ReferenceName());
}

TEST(CompareTest, ReferenceStartOk)
{
    BamRecord r1; r1.Impl().Position(30);
    BamRecord r2; r2.Impl().Position(20);
    BamRecord r3; r3.Impl().Position(40);
    BamRecord r4; r4.Impl().Position(10);

    auto records = std::vector<BamRecord>{ r1, r2, r3, r4 };
    std::sort(records.begin(), records.end(), Compare::ReferenceStart());

    EXPECT_EQ(r4.ReferenceStart(), records.at(0).ReferenceStart());
    EXPECT_EQ(r2.ReferenceStart(), records.at(1).ReferenceStart());
    EXPECT_EQ(r1.ReferenceStart(), records.at(2).ReferenceStart());
    EXPECT_EQ(r3.ReferenceStart(), records.at(3).ReferenceStart());
}

TEST(CompareTest, ZmwOk)
{
    int32_t z1 = 30;
    int32_t z2 = 20;
    int32_t z3 = 40;
    int32_t z4 = 10;

    auto records = std::vector<BamRecord>
    {
        CompareTests::makeRecordWithTag("zm", Tag(z1)),
        CompareTests::makeRecordWithTag("zm", Tag(z2)),
        CompareTests::makeRecordWithTag("zm", Tag(z3)),
        CompareTests::makeRecordWithTag("zm", Tag(z4))
    };
    std::sort(records.begin(), records.end(), Compare::Zmw());

    EXPECT_EQ(z4, records.at(0).HoleNumber());
    EXPECT_EQ(z2, records.at(1).HoleNumber());
    EXPECT_EQ(z1, records.at(2).HoleNumber());
    EXPECT_EQ(z3, records.at(3).HoleNumber());
}

// clang-format on

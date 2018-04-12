// Author: Derek Barnett

#include <chrono>
#include <cstdint>
#include <string>

#include <gtest/gtest.h>

#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordView.h>
#include <pbbam/BamTagCodec.h>

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

typedef vector<uint16_t> f_data;

namespace BamRecordClippingTests {

static 
ReadGroupInfo MakeReadGroup(const FrameCodec codec,
                            const std::string& movieName,
                            const std::string& readType)
{
    ReadGroupInfo rg{movieName, readType};
    rg.IpdCodec(codec);
    rg.PulseWidthCodec(codec);
    return rg;
}

static
BamRecord MakeRecord(const Position qStart,
                     const Position qEnd,
                     const string& seq,
                     const string& quals,
                     const string& tagBases,
                     const string& tagQuals,
                     const f_data& frames,
                     const string& pulseCall = "",
                     const string& pulseBases = "",
                     const string& pulseQuals = "",
                     const f_data& pulseFrames = f_data(),
                     const FrameCodec codec = FrameCodec::RAW)
{
    BamRecordImpl impl;
    impl.SetSequenceAndQualities(seq, quals);

    TagCollection tags;
    tags["qs"] = qStart;        // qStart
    tags["qe"] = qEnd;          // qEnd
    tags["dt"] = tagBases;      // deletionTag
    tags["st"] = tagBases;      // substitutionTag
    tags["dq"] = tagQuals;      // deletionQV
    tags["iq"] = tagQuals;      // insertionQV
    tags["mq"] = tagQuals;      // mergeQV
    tags["sq"] = tagQuals;      // substitutionQV
    tags["ip"] = frames;        // IPD
    tags["pw"] = frames;        // pulseWidth
    tags["pc"] = pulseCall;     // pulseCall
    tags["pt"] = pulseBases;    // altLabelTag
    tags["pq"] = pulseQuals;    // labelQV
    tags["pv"] = pulseQuals;    // altLabelQV
    tags["pg"] = pulseQuals;    // pulseMergeQV
    tags["pa"] = pulseFrames;   // pkmean
    tags["pm"] = pulseFrames;   // pkmid
    impl.Tags(tags);

    const auto rg = MakeReadGroup(codec, "movie", "SUBREAD");

    BamRecord bam(std::move(impl));
    bam.header_.AddReadGroup(rg);
    bam.ReadGroup(rg);
    return bam;
}

static
BamRecord MakeCCSRecord(const string& seq,
                        const string& quals,
                        const string& tagBases,
                        const string& tagQuals,
                        const f_data& frames,
                        const string& pulseCall = "",
                        const string& pulseBases = "",
                        const string& pulseQuals = "",
                        const f_data& pulseFrames = f_data(),
                        const FrameCodec codec = FrameCodec::RAW)
{
    BamRecordImpl impl;
    impl.Name("movie/42/ccs");
    impl.SetSequenceAndQualities(seq, quals);

    TagCollection tags;
    tags["dt"] = tagBases;      // deletionTag
    tags["st"] = tagBases;      // substitutionTag
    tags["dq"] = tagQuals;      // deletionQV
    tags["iq"] = tagQuals;      // insertionQV
    tags["mq"] = tagQuals;      // mergeQV
    tags["sq"] = tagQuals;      // substitutionQV
    tags["ip"] = frames;        // IPD
    tags["pw"] = frames;        // pulseWidth
    tags["pc"] = pulseCall;     // pulseCall
    tags["pt"] = pulseBases;    // altLabelTag
    tags["pq"] = pulseQuals;    // labelQV
    tags["pv"] = pulseQuals;    // altLabelQV
    tags["pg"] = pulseQuals;    // pulseMergeQV
    tags["pa"] = pulseFrames;   // pkmean
    tags["pm"] = pulseFrames;   // pkmid
    impl.Tags(tags);

    const auto rg = MakeReadGroup(codec, "movie", "CCS");

    BamRecord bam(std::move(impl));
    bam.header_.AddReadGroup(rg);
    bam.ReadGroup(rg);
    return bam;
}

} // namespace BamRecordClippingTests

TEST(BamRecordClippingTest, ClipToQuery_Basic)
{
    const Position qStart  = 500;
    const Position qEnd    = 510;
    const string seq       = "AACCGTTAGC";
    const string quals     = "?]?]?]?]?*";
    const string tagBases  = "AACCGTTAGC";
    const string tagQuals  = "?]?]?]?]?*";
    const f_data frames    = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const string pulseCall   = "ttAaAtaCCGggatTTAcatGCt";
    const string pulseBases  = pulseCall;
    const string pulseQuals  = "==?=]==?]?====]?]===?*=";
    const f_data pulseFrames = { 0,0,10,0,10,0,0,20,20,30,0,0,0,0,40,40,10,0,0,0,30,20,0 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const string seq_clipped      = "CCGTTAG";
    const string quals_clipped    = "?]?]?]?";
    const string tagBases_clipped = "CCGTTAG";
    const string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const string pulseCall_clipped = "CCGggatTTAcatG";
    const string pulseQuals_clipped = "?]?====]?]===?";
    const f_data pulseFrames_clipped = { 20,20,30,0,0,0,0,40,40,10,0,0,0,30 };

    const string seq_rev       = "GCTAACGGTT";
    const string pulseCall_rev = "aGCatgTAAatccCGGtaTtTaa";
    const string quals_rev     = "*?]?]?]?]?";
    const string tagQuals_rev  = quals_rev;
    const f_data frames_rev    = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const string seq_rev_clipped   = "CTAACGG";
    const string quals_rev_clipped = "?]?]?]?";
    const string tagBases_rev_clipped = seq_rev_clipped;
    const string tagQuals_rev_clipped = quals_rev_clipped;
    const f_data frames_rev_clipped = { 30, 10, 40, 40, 30, 20, 20 };

    const string pulseCall_rev_clipped = "CatgTAAatccCGG";
    const string pulseQuals_rev_clipped    = "?===]?]====?]?";
    const f_data pulseFrames_rev_clipped = { 30,0,0,0,10,40,40,0,0,0,0,30,20,20 };

    const string s1_cigar = "10=";
    const string s2_cigar = "5=3D5=";
    const string s3_cigar = "4=1D2I2D4=";

    const string s1_cigar_clipped = "7=";
    const string s2_cigar_clipped = "3=3D4=";
    const string s3_cigar_clipped = "2=1D2I2D3=";

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  pulseCall, pulseBases, pulseQuals, pulseFrames);

    BamRecord s0 = prototype; // unmapped record
    BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
    BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
    BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
    BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
    BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
    BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);

    s0.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s1.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    {   // s0

        EXPECT_FALSE(s0.IsMapped());
        EXPECT_EQ(clipStart, s0.QueryStart());
        EXPECT_EQ(clipEnd,   s0.QueryEnd());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.AlignedStart());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.AlignedEnd());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.ReferenceStart());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.ReferenceEnd());

        const BamRecordView view
        {
            s0,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,       view.Sequence());
        EXPECT_EQ(quals_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_clipped, view.PulseCalls());
    }

    {   // s1 - FORWARD

        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(clipStart, s1.QueryStart());
        EXPECT_EQ(clipEnd,   s1.QueryEnd());
        EXPECT_EQ(clipStart, s1.AlignedStart());   // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s1.AlignedEnd());     // alignStart + seqLength
        EXPECT_EQ(102, s1.ReferenceStart());       // 100 + startOffset
        EXPECT_EQ(109, s1.ReferenceEnd());         // RefStart + 7=

        EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

        const BamRecordView view
        {
            s1,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,       view.Sequence());
        EXPECT_EQ(quals_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_clipped, view.PulseCalls());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s1_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s1_rev.QueryEnd());
        EXPECT_EQ(clipStart, s1_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s1_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(102, s1_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(109, s1_rev.ReferenceEnd());          // RefStart + 7=

        EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_rev_clipped,       view.Sequence());
        EXPECT_EQ(quals_rev_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_rev_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_rev_clipped, view.PulseCalls());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(clipStart, s2.QueryStart());
        EXPECT_EQ(clipEnd,   s2.QueryEnd());
        EXPECT_EQ(clipStart, s2.AlignedStart());   // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s2.AlignedEnd());     // alignStart + seqLength
        EXPECT_EQ(102, s2.ReferenceStart());       // 100 + startOffset
        EXPECT_EQ(112, s2.ReferenceEnd());         // RefStart + 7= + 3D

        EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

        const BamRecordView view
        {
            s2,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,      view.Sequence());
        EXPECT_EQ(quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s2_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
        EXPECT_EQ(clipStart, s2_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s2_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(102, s2_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(112, s2_rev.ReferenceEnd());          // RefStart + 7= + 3D

        EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_rev_clipped,       view.Sequence());
        EXPECT_EQ(quals_rev_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_rev_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_rev_clipped, view.PulseCalls());
    }

    {   // s3 - FORWARD

        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(clipStart, s3.QueryStart());
        EXPECT_EQ(clipEnd,   s3.QueryEnd());
        EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
        EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5= + 3D

        EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

        const BamRecordView view
        {
            s3,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,       view.Sequence());
        EXPECT_EQ(quals_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_clipped, view.PulseCalls());
    }

    {   // s3 - REVERSE

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s3_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
        EXPECT_EQ(clipStart, s3_rev.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s3_rev.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(102, s3_rev.ReferenceStart());         // 100 + startOffset
        EXPECT_EQ(110, s3_rev.ReferenceEnd());           // RefStart + 5= + 3D

        EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_rev_clipped,       view.Sequence());
        EXPECT_EQ(quals_rev_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_rev_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_rev_clipped, view.PulseCalls());
    }
}

TEST(BamRecordClippingTest, ClipToQuery_WithSoftClips)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const string seq      = "TTAACCGTTAGCAAA";
    const string seq_rev  = "TTTGCTAACGGTTAA";
    const string quals    = "--?]?]?]?]?*+++";
    const string tagBases = "TTAACCGTTAGCAAA";
    const string tagQuals = "--?]?]?]?]?*+++";
    const string tagQuals_rev = "+++*?]?]?]?]?--";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };
    const f_data frames_rev = { 10, 10, 10, 20, 30, 10, 40, 40, 30, 20, 20, 10, 10, 40, 40 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const string s1_cigar = "2S10=3S";
    const string s1_cigar_clipped = "7=";
    const string s1_seq_clipped      = "AACCGTT";
    const string s1_quals_clipped    = "?]?]?]?";
    const string s1_tagBases_clipped = s1_seq_clipped;
    const string s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const string s1_seq_rev_clipped   = "AACGGTT";
    const string s1_quals_rev_clipped = "?]?]?]?";
    const string s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const string s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 40, 40, 30, 20, 20, 10, 10 };

    const string s2_cigar = "2S5=3D5=3S";
    const string s2_cigar_clipped = "5=3D2=";
    const string s2_seq_clipped      = "AACCGTT";
    const string s2_quals_clipped    = "?]?]?]?";
    const string s2_tagBases_clipped = s2_seq_clipped;
    const string s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const string s2_seq_rev_clipped   = "AACGGTT";
    const string s2_quals_rev_clipped = "?]?]?]?";
    const string s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const string s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 40, 40, 30, 20, 20, 10, 10 };

    const string s3_cigar = "2S4=1D2I2D4=3S";
    const string s3_cigar_clipped = "4=1D2I2D1=";
    const string s3_seq_clipped      = "AACCGTT";
    const string s3_quals_clipped    = "?]?]?]?";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const string s3_seq_rev_clipped   = "AACGGTT";
    const string s3_quals_rev_clipped = "?]?]?]?";
    const string s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const string s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
    const f_data s3_frames_rev_clipped = { 40, 40, 30, 20, 20, 10, 10 };

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  seq, tagBases, tagQuals, frames);
    BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
    BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
    BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
    BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
    BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
    BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);

    // sanity checks before clipping
    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(tPos, s1.ReferenceStart());
    EXPECT_EQ(tPos + 10, s1.ReferenceEnd()); // 10=

    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(tPos, s1_rev.ReferenceStart());
    EXPECT_EQ(tPos + 10, s1_rev.ReferenceEnd()); // 10=

    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(tPos, s2.ReferenceStart());
    EXPECT_EQ(tPos + 13, s2.ReferenceEnd());   // 5= + 3D + 5=

    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(tPos, s2_rev.ReferenceStart());
    EXPECT_EQ(tPos + 13, s2_rev.ReferenceEnd());   // 5= + 3D + 5=

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(tPos, s3.ReferenceStart());
    EXPECT_EQ(tPos + 11, s3.ReferenceEnd());   // 4= + 1D + 2D + 4=

    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(tPos, s3_rev.ReferenceStart());
    EXPECT_EQ(tPos + 11, s3_rev.ReferenceEnd());   // 4= + 1D + 2D + 4=

    s1.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    {   // s1 - FORWARD

        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(clipStart, s1.QueryStart());
        EXPECT_EQ(clipEnd,   s1.QueryEnd());
        EXPECT_EQ(clipStart, s1.AlignedStart());    // queryStart (no soft clips left)
        EXPECT_EQ(clipEnd,   s1.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s1.ReferenceStart());  // tPos
        EXPECT_EQ(tPos + 7,  s1.ReferenceEnd());    // RefStart + 7=

        EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

        const BamRecordView view
        {
            s1,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s1_seq_clipped,      view.Sequence());
        EXPECT_EQ(s1_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s1_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s1_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s1_frames_clipped,   view.IPD().Data());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s1_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s1_rev.QueryEnd());
        EXPECT_EQ(clipStart, s1_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s1_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s1_rev.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 7,  s1_rev.ReferenceEnd());    // RefStart + 7=

        EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s1_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s1_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s1_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s1_frames_rev_clipped,   view.IPD().Data());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(clipStart, s2.QueryStart());
        EXPECT_EQ(clipEnd,   s2.QueryEnd());
        EXPECT_EQ(clipStart, s2.AlignedStart());    // queryStart (no soft clips left)
        EXPECT_EQ(clipEnd,   s2.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s2.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 10, s2.ReferenceEnd());    // RefStart + 5=3D2=

        EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

        const BamRecordView view
        {
            s2,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s2_seq_clipped,      view.Sequence());
        EXPECT_EQ(s2_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s2_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s2_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s2_frames_clipped,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s2_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
        EXPECT_EQ(clipStart, s2_rev.AlignedStart());    // queryStart (no soft clips left)
        EXPECT_EQ(clipEnd,   s2_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s2_rev.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 10, s2_rev.ReferenceEnd());    // RefStart + 5=3D2=

        EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s2_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s2_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s2_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s2_frames_rev_clipped,   view.IPD().Data());
    }

    {   // s3 - FORWARD

        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(clipStart, s3.QueryStart());
        EXPECT_EQ(clipEnd,   s3.QueryEnd());
        EXPECT_EQ(clipStart, s3.AlignedStart());    // queryStart (no soft clips left)
        EXPECT_EQ(clipEnd,   s3.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s3.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 8,  s3.ReferenceEnd());    // RefStart + 4=1D2D1=

        EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

        const BamRecordView view
        {
            s3,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s3_seq_clipped,      view.Sequence());
        EXPECT_EQ(s3_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s3_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s3_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s3_frames_clipped,   view.IPD().Data());
    }

    {   // s3 - REVERSE
        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s3_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
        EXPECT_EQ(clipStart, s3_rev.AlignedStart());    // queryStart (no soft clips left)
        EXPECT_EQ(clipEnd,   s3_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s3_rev.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 8,  s3_rev.ReferenceEnd());    // RefStart + 4=1D2D1=

        EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s3_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s3_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s3_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s3_frames_rev_clipped,   view.IPD().Data());
    }
}

TEST(BamRecordClippingTest, ClipToReference_Basic)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const string tagQuals_rev = "*?]?]?]?]?";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const string s1_cigar = "10=";
    const string s1_cigar_clipped = "5=";
    const string s1_seq_clipped      = "CCGTT";
    const string s1_quals_clipped    = "?]?]?";
    const string s1_tagBases_clipped = s1_seq_clipped;
    const string s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 20, 20, 30, 40, 40 };
    const string s1_seq_rev_clipped   = "TAACG";
    const string s1_quals_rev_clipped = "]?]?]";
    const string s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const string s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 10, 40, 40, 30, 20 };

    const string s2_cigar = "5=3D5=";
    const string s2_cigar_clipped = "3=2D";
    const string s2_seq_clipped      = "CCG";
    const string s2_quals_clipped    = "?]?";
    const string s2_tagBases_clipped = s2_seq_clipped;
    const string s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 20, 20, 30 };
    const string s2_seq_rev_clipped   = "TAA";
    const string s2_quals_rev_clipped = "]?]";
    const string s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const string s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 10, 40, 40 };

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };
    const string s3_seq_rev_clipped   = "TAAC";
    const string s3_quals_rev_clipped = "]?]?";
    const string s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const string s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
    const f_data s3_frames_rev_clipped = { 10, 40, 40, 30};

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  seq, tagBases, tagQuals, frames);
    BamRecord s0 = prototype;
    BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
    BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
    BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
    BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
    BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
    BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);

    s0.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s1.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);

    {   // s0 - no clipping should have been done to unmapped record

        EXPECT_FALSE(s0.IsMapped());
        EXPECT_EQ(prototype.QueryStart(),     s0.QueryStart());
        EXPECT_EQ(prototype.QueryEnd(),       s0.QueryEnd());
        EXPECT_EQ(prototype.AlignedStart(),   s0.AlignedStart());
        EXPECT_EQ(prototype.AlignedEnd(),     s0.AlignedEnd());
        EXPECT_EQ(prototype.ReferenceStart(), s0.ReferenceStart());
        EXPECT_EQ(prototype.ReferenceEnd(),   s0.ReferenceEnd());

        const BamRecordView protoView
        {
            prototype,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        const BamRecordView view
        {
            s0,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(protoView.Sequence(),       view.Sequence());
        EXPECT_EQ(protoView.Qualities(),      view.Qualities());
        EXPECT_EQ(protoView.DeletionTags(),    view.DeletionTags());
        EXPECT_EQ(protoView.DeletionQVs(),     view.DeletionQVs());
        EXPECT_EQ(protoView.LabelQVs(),        view.LabelQVs());
        EXPECT_EQ(protoView.AltLabelQVs(),     view.AltLabelQVs());
        EXPECT_EQ(protoView.IPD(),            view.IPD());
    }

    {   // s1 - FORWARD

        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(502,   s1.QueryStart());
        EXPECT_EQ(507,   s1.QueryEnd());
        EXPECT_EQ(502,   s1.AlignedStart());       // queryStart (no soft clips)
        EXPECT_EQ(507,   s1.AlignedEnd());         // alignStart + seqLength
        EXPECT_EQ(clipStart, s1.ReferenceStart()); // clipStart
        EXPECT_EQ(clipEnd,   s1.ReferenceEnd());   // clipEnd

        EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

        const BamRecordView view
        {
            s1,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s1_seq_clipped,      view.Sequence());
        EXPECT_EQ(s1_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s1_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s1_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s1_frames_clipped,   view.IPD().Data());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(503, s1_rev.QueryStart());
        EXPECT_EQ(508, s1_rev.QueryEnd());
        EXPECT_EQ(503, s1_rev.AlignedStart());          // queryStart (no soft clips)
        EXPECT_EQ(508, s1_rev.AlignedEnd());            // alignStart + seqLength
        EXPECT_EQ(clipStart, s1_rev.ReferenceStart());  // clipStart
        EXPECT_EQ(clipEnd,   s1_rev.ReferenceEnd());    // clipEnd

        EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s1_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s1_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s1_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s1_frames_rev_clipped,   view.IPD().Data());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(502, s2.QueryStart());
        EXPECT_EQ(505, s2.QueryEnd());
        EXPECT_EQ(502, s2.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(505, s2.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(clipStart, s2.ReferenceStart());   // clipStart
        EXPECT_EQ(clipEnd,   s2.ReferenceEnd());     // clipEnd

        EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

        const BamRecordView view
        {
            s2,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s2_seq_clipped,      view.Sequence());
        EXPECT_EQ(s2_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s2_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s2_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s2_frames_clipped,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(505, s2_rev.QueryStart());
        EXPECT_EQ(508, s2_rev.QueryEnd());
        EXPECT_EQ(505, s2_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(508, s2_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(clipStart, s2_rev.ReferenceStart());  // clipStart
        EXPECT_EQ(clipEnd,   s2_rev.ReferenceEnd());    // clipEnd

        EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s2_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s2_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s2_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s2_frames_rev_clipped,   view.IPD().Data());
    }

    {   // s3 - FORWARD

        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(502, s3.QueryStart());
        EXPECT_EQ(506, s3.QueryEnd());
        EXPECT_EQ(502, s3.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(506, s3.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(clipStart, s3.ReferenceStart());   // clipStart
        EXPECT_EQ(clipEnd,   s3.ReferenceEnd());     // clipEnd

        EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

        const BamRecordView view
        {
            s3,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s3_seq_clipped,      view.Sequence());
        EXPECT_EQ(s3_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s3_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s3_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s3_frames_clipped,   view.IPD().Data());
    }

    {   // s3 - REVERSE

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(504, s3_rev.QueryStart());
        EXPECT_EQ(508, s3_rev.QueryEnd());
        EXPECT_EQ(504, s3_rev.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(508, s3_rev.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(clipStart, s3_rev.ReferenceStart());  // clipStart
        EXPECT_EQ(clipEnd,   s3_rev.ReferenceEnd());    // clipEnd

        EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s3_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s3_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s3_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s3_frames_rev_clipped,   view.IPD().Data());
    }
}

TEST(BamRecordClippingTest, ClipToReference_WithSoftClips)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const string seq      = "TTAACCGTTAGCAAA";
    const string quals    = "--?]?]?]?]?*+++";
    const string tagBases = "TTAACCGTTAGCAAA";
    const string tagQuals = "--?]?]?]?]?*+++";
    const string tagQuals_rev = "+++*?]?]?]?]?--";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const string seq_rev      = "TTTGCTAACGGTTAA";
    const string quals_rev    = "+++*?]?]?]?]?--";
    const f_data frames_rev   = { 10, 10, 10, 20, 30, 10, 40, 40, 30, 20, 20, 10, 10, 40, 40 };

    const string s1_cigar = "2S10=3S";
    const string s1_cigar_clipped = "5=";
    const string s1_seq_clipped      = "CCGTT";
    const string s1_quals_clipped    = "?]?]?";
    const string s1_tagBases_clipped = s1_seq_clipped;
    const string s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 20, 20, 30, 40, 40 };
    const string s1_seq_rev_clipped   = "CTAAC";
    const string s1_quals_rev_clipped = "?]?]?";
    const string s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const string s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 30, 10, 40, 40, 30 };

    const string s2_cigar = "2S5=3D5=3S";
    const string s2_cigar_clipped = "3=2D";
    const string s2_seq_clipped      = "CCG";
    const string s2_quals_clipped    = "?]?";
    const string s2_tagBases_clipped = s2_seq_clipped;
    const string s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 20, 20, 30 };
    const string s2_seq_rev_clipped   = "CTA";
    const string s2_quals_rev_clipped = "?]?";
    const string s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const string s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 30, 10, 40 };

    const string s3_cigar = "2S4=1D2I2D4=3S";
    const string s3_cigar_clipped = "2=1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };
    const string s3_seq_rev_clipped   = "CTAA";
    const string s3_quals_rev_clipped = "?]?]";
    const string s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const string s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
    const f_data s3_frames_rev_clipped = { 30, 10, 40, 40 };

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  seq, tagBases, tagQuals, frames);
    BamRecord s0 = prototype;
    BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
    BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
    BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
    BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
    BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
    BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);

    // sanity checks before clipping
    EXPECT_FALSE(s0.IsMapped());

    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(500,       s1.QueryStart());      // queryStart
    EXPECT_EQ(515,       s1.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(502,       s1.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(512,       s1.AlignedEnd());      // alignedStart + 10=
    EXPECT_EQ(tPos,      s1.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 10, s1.ReferenceEnd());    // tPos + 10=

    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(500,       s1_rev.QueryStart());      // queryStart
    EXPECT_EQ(515,       s1_rev.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(503,       s1_rev.AlignedStart());    // queryStart + 3S
    EXPECT_EQ(513,       s1_rev.AlignedEnd());      // alignedStart + 10=
    EXPECT_EQ(tPos,      s1_rev.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 10, s1_rev.ReferenceEnd());    // tPos + 10=

    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(500,       s2.QueryStart());      // queryStart
    EXPECT_EQ(515,       s2.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(502,       s2.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(512,       s2.AlignedEnd());      // alignedStart + 5=5=
    EXPECT_EQ(tPos,      s2.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 13, s2.ReferenceEnd());    // tPos + 5=3D5=

    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(500,       s2_rev.QueryStart());      // queryStart
    EXPECT_EQ(515,       s2_rev.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(503,       s2_rev.AlignedStart());    // queryStart + S
    EXPECT_EQ(513,       s2_rev.AlignedEnd());      // alignedStart + 5=5=
    EXPECT_EQ(tPos,      s2_rev.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 13, s2_rev.ReferenceEnd());    // tPos + 5=3D5=

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(500,       s3.QueryStart());      // queryStart
    EXPECT_EQ(515,       s3.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(502,       s3.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(512,       s3.AlignedEnd());      // alignedStart + 4=2I4=
    EXPECT_EQ(tPos,      s3.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 11, s3.ReferenceEnd());    // tPos + 4=1D2D4=

    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(500,       s3_rev.QueryStart());      // queryStart
    EXPECT_EQ(515,       s3_rev.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(503,       s3_rev.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(513,       s3_rev.AlignedEnd());      // alignedStart + 4=2I4=
    EXPECT_EQ(tPos,      s3_rev.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 11, s3_rev.ReferenceEnd());    // tPos + 4=1D2D4=

    s0.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s1.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);

    {   // s0 - no clipping should have been done to unmapped record

        EXPECT_FALSE(s0.IsMapped());
        EXPECT_EQ(prototype.QueryStart(),     s0.QueryStart());
        EXPECT_EQ(prototype.QueryEnd(),       s0.QueryEnd());
        EXPECT_EQ(prototype.AlignedStart(),   s0.AlignedStart());
        EXPECT_EQ(prototype.AlignedEnd(),     s0.AlignedEnd());
        EXPECT_EQ(prototype.ReferenceStart(), s0.ReferenceStart());
        EXPECT_EQ(prototype.ReferenceEnd(),   s0.ReferenceEnd());

        const BamRecordView protoView
        {
            prototype,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        const BamRecordView view
        {
            s0,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(protoView.Sequence(),      view.Sequence());
        EXPECT_EQ(protoView.Qualities(),     view.Qualities());
        EXPECT_EQ(protoView.DeletionTags(),  view.DeletionTags());
        EXPECT_EQ(protoView.DeletionQVs(),   view.DeletionQVs());
        EXPECT_EQ(protoView.LabelQVs(),      view.LabelQVs());
        EXPECT_EQ(protoView.AltLabelQVs(),   view.AltLabelQVs());
        EXPECT_EQ(protoView.IPD(),           view.IPD());
    }

    {   // s1 - FORWARD

        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(504,   s1.QueryStart());         // new queryStart
        EXPECT_EQ(509,   s1.QueryEnd());           // queryStart + new seqLength
        EXPECT_EQ(504,   s1.AlignedStart());       // queryStart (no soft clips remaining)
        EXPECT_EQ(509,   s1.AlignedEnd());         // alignStart + new seqLength
        EXPECT_EQ(clipStart, s1.ReferenceStart()); // clipStart
        EXPECT_EQ(clipEnd,   s1.ReferenceEnd());   // clipEnd

        EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

        const BamRecordView view
        {
            s1,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s1_seq_clipped,      view.Sequence());
        EXPECT_EQ(s1_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s1_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s1_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s1_frames_clipped,   view.IPD().Data());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(506,   s1_rev.QueryStart());         // new queryStart
        EXPECT_EQ(511,   s1_rev.QueryEnd());           // queryStart + new seqLength
        EXPECT_EQ(506,   s1_rev.AlignedStart());       // queryStart (no soft clips remaining)
        EXPECT_EQ(511,   s1_rev.AlignedEnd());         // alignStart + new seqLength
        EXPECT_EQ(clipStart, s1_rev.ReferenceStart()); // clipStart
        EXPECT_EQ(clipEnd,   s1_rev.ReferenceEnd());   // clipEnd

        EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s1_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s1_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s1_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s1_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s1_frames_rev_clipped,   view.IPD().Data());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(504, s2.QueryStart());
        EXPECT_EQ(507, s2.QueryEnd());
        EXPECT_EQ(504, s2.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(507, s2.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(clipStart, s2.ReferenceStart());   // clipStart
        EXPECT_EQ(clipEnd,   s2.ReferenceEnd());     // clipEnd

        EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

        const BamRecordView view
        {
            s2,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s2_seq_clipped,      view.Sequence());
        EXPECT_EQ(s2_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s2_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s2_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s2_frames_clipped,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(508,   s2_rev.QueryStart());         // new queryStart
        EXPECT_EQ(511,   s2_rev.QueryEnd());           // queryStart + new seqLength
        EXPECT_EQ(508,   s2_rev.AlignedStart());       // queryStart (no soft clips remaining)
        EXPECT_EQ(511,   s2_rev.AlignedEnd());         // alignStart + new seqLength
        EXPECT_EQ(clipStart, s2_rev.ReferenceStart()); // clipStart
        EXPECT_EQ(clipEnd,   s2_rev.ReferenceEnd());   // clipEnd

        EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s2_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s2_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s2_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s2_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s2_frames_rev_clipped,   view.IPD().Data());
    }

    {   // s3 - FORWARD
        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(504, s3.QueryStart());
        EXPECT_EQ(508, s3.QueryEnd());
        EXPECT_EQ(504, s3.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(508, s3.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(clipStart, s3.ReferenceStart());   // clipStart
        EXPECT_EQ(clipEnd,   s3.ReferenceEnd());     // clipEnd

        EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

        const BamRecordView view
        {
            s3,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s3_seq_clipped,      view.Sequence());
        EXPECT_EQ(s3_quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s3_tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(s3_tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s3_frames_clipped,   view.IPD().Data());
    }

    {   // s3 - REVERSE
        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(507,   s3_rev.QueryStart());         // new queryStart
        EXPECT_EQ(511,   s3_rev.QueryEnd());           // queryStart + new seqLength
        EXPECT_EQ(507,   s3_rev.AlignedStart());       // queryStart (no soft clips remaining)
        EXPECT_EQ(511,   s3_rev.AlignedEnd());         // alignStart + new seqLength
        EXPECT_EQ(clipStart, s3_rev.ReferenceStart()); // clipStart
        EXPECT_EQ(clipEnd,   s3_rev.ReferenceEnd());   // clipEnd

        EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(s3_seq_rev_clipped,      view.Sequence());
        EXPECT_EQ(s3_quals_rev_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(s3_tagBases_rev_clipped, view.DeletionTags());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(s3_tagQuals_rev_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(s3_frames_rev_clipped,   view.IPD().Data());
    }
}

TEST(BamRecordClippingTest, ClippedToQueryCopy)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const string seq_clipped      = "CCGTTAG";
    const string quals_clipped    = "?]?]?]?";
    const string tagBases_clipped = "CCGTTAG";
    const string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D3=";

    BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                            seq, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    BamRecord s3 = prototype.Clipped(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(clipStart, s3.QueryStart());
    EXPECT_EQ(clipEnd,   s3.QueryEnd());
    EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5= + 3D

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    const BamRecordView view
    {
        s3,
        Orientation::GENOMIC,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(seq_clipped,      view.Sequence());
    EXPECT_EQ(quals_clipped,    view.Qualities().Fastq());
    EXPECT_EQ(tagBases_clipped, view.DeletionTags());
    EXPECT_EQ(tagQuals_clipped, view.DeletionQVs().Fastq());
    EXPECT_EQ(tagQuals_clipped, view.LabelQVs().Fastq());
    EXPECT_EQ(tagQuals_clipped, view.AltLabelQVs().Fastq());
    EXPECT_EQ(frames_clipped,   view.IPD().Data());
}

TEST(BamRecordClippingTest, ClippedToReferenceCopy)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };

    BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                            seq, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    const BamRecord s3 = prototype.Clipped(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);

    // s3 - FORWARD
    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(502, s3.QueryStart());
    EXPECT_EQ(506, s3.QueryEnd());
    EXPECT_EQ(502, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(506, s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(clipStart, s3.ReferenceStart());   // clipStart
    EXPECT_EQ(clipEnd,   s3.ReferenceEnd());     // clipEnd

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    const BamRecordView view
    {
        s3,
        Orientation::GENOMIC,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(s3_seq_clipped,      view.Sequence());
    EXPECT_EQ(s3_quals_clipped,    view.Qualities().Fastq());
    EXPECT_EQ(s3_tagBases_clipped, view.DeletionTags());
    EXPECT_EQ(s3_tagQuals_clipped, view.DeletionQVs().Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, view.LabelQVs().Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, view.AltLabelQVs().Fastq());
    EXPECT_EQ(s3_frames_clipped,   view.IPD().Data());
}

TEST(BamRecordClippingTest, StaticClippedToQuery)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const string seq_clipped      = "CCGTTAG";
    const string quals_clipped    = "?]?]?]?";
    const string tagBases_clipped = "CCGTTAG";
    const string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D3=";

    BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                            seq, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    const BamRecord s3 = BamRecord::Clipped(prototype, ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(clipStart, s3.QueryStart());
    EXPECT_EQ(clipEnd,   s3.QueryEnd());
    EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5= + 3D

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    const BamRecordView view
    {
        s3,
        Orientation::GENOMIC,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(seq_clipped,      view.Sequence());
    EXPECT_EQ(quals_clipped,    view.Qualities().Fastq());
    EXPECT_EQ(tagBases_clipped, view.DeletionTags());
    EXPECT_EQ(tagQuals_clipped, view.DeletionQVs().Fastq());
    EXPECT_EQ(tagQuals_clipped, view.LabelQVs().Fastq());
    EXPECT_EQ(tagQuals_clipped, view.AltLabelQVs().Fastq());
    EXPECT_EQ(frames_clipped,   view.IPD().Data());
}

TEST(BamRecordClippingTest, StaticClippedToReference)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };

    BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                            seq, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    const BamRecord s3 = BamRecord::Clipped(prototype, ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);

    // s3 - FORWARD
    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(502, s3.QueryStart());
    EXPECT_EQ(506, s3.QueryEnd());
    EXPECT_EQ(502, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(506, s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(clipStart, s3.ReferenceStart());   // clipStart
    EXPECT_EQ(clipEnd,   s3.ReferenceEnd());     // clipEnd

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    const BamRecordView view
    {
        s3,
        Orientation::GENOMIC,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(s3_seq_clipped,      view.Sequence());
    EXPECT_EQ(s3_quals_clipped,    view.Qualities().Fastq());
    EXPECT_EQ(s3_tagBases_clipped, view.DeletionTags());
    EXPECT_EQ(s3_tagQuals_clipped, view.DeletionQVs().Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, view.LabelQVs().Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, view.AltLabelQVs().Fastq());
    EXPECT_EQ(s3_frames_clipped,   view.IPD().Data());
}

TEST(BamRecordTest, ClipCigarData)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const string seq      = "TTAACCGTTAGCAAA";
    const string quals    = "--?]?]?]?]?*+++";
    const string tagBases = "TTAACCGTTAGCAAA";
    const string tagQuals = "--?]?]?]?]?*+++";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };
    const uint8_t mapQual = 80;
    BamRecord s3 = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                     seq, tagBases, tagQuals, frames);
    BamRecord s3_rev = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                         seq, tagBases, tagQuals, frames);

    const string s3_cigar = "5H2S4=1D2I2D4=3S7H";
    s3.Map(0, 100, Strand::FORWARD, s3_cigar, mapQual);
    s3_rev.Map(0, 100, Strand::REVERSE, s3_cigar, mapQual);

    const Cigar s3_cigar_raw     = s3.CigarData();
    const Cigar s3_cigar_clipped = s3.CigarData(true);

    EXPECT_EQ(s3_cigar, s3_cigar_raw.ToStdString());
    EXPECT_EQ(string("4=1D2I2D4="), s3_cigar_clipped.ToStdString());
}

TEST(BamRecordTest, CCS_ClipToQuery)
{
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 2;
    const Position clipEnd   = 9;

    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const string seq_clipped      = "CCGTTAG";
    const string quals_clipped    = "?]?]?]?";
    const string tagBases_clipped = "CCGTTAG";
    const string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D3=";

    BamRecord prototype = BamRecordClippingTests::MakeCCSRecord(seq, quals, tagBases, tagQuals, frames,
                                               seq, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    BamRecord s3 = prototype.Clipped(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(0,   s3.AlignedStart());     // record start (no soft clips)
    EXPECT_EQ(7,   s3.AlignedEnd());       // alignStart + clipped seqLength
    EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5= + 3D

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    const BamRecordView view
    {
        s3,
        Orientation::GENOMIC,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(seq_clipped,      view.Sequence());
    EXPECT_EQ(quals_clipped,    view.Qualities().Fastq());
    EXPECT_EQ(tagBases_clipped, view.DeletionTags());
    EXPECT_EQ(tagQuals_clipped, view.DeletionQVs().Fastq());
    EXPECT_EQ(tagQuals_clipped, view.LabelQVs().Fastq());
    EXPECT_EQ(tagQuals_clipped, view.AltLabelQVs().Fastq());
    EXPECT_EQ(frames_clipped,   view.IPD().Data());
}

TEST(BamRecordTest, CCS_ClipToReference)
{
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const string s3_cigar = "4=1D2I2D4=";
    const string s3_cigar_clipped = "2=1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };

    BamRecord prototype = BamRecordClippingTests::MakeCCSRecord(seq, quals, tagBases, tagQuals, frames,
                                               seq, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    const BamRecord s3 = BamRecord::Clipped(prototype, ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(0, s3.AlignedStart());     // record tart (no soft clips)
    EXPECT_EQ(4, s3.AlignedEnd());       // alignStart + clipped seqLength (4)
    EXPECT_EQ(clipStart, s3.ReferenceStart());   // clipStart
    EXPECT_EQ(clipEnd,   s3.ReferenceEnd());     // clipEnd

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    const BamRecordView view
    {
        s3,
        Orientation::GENOMIC,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(s3_seq_clipped,      view.Sequence());
    EXPECT_EQ(s3_quals_clipped,    view.Qualities().Fastq());
    EXPECT_EQ(s3_tagBases_clipped, view.DeletionTags());
    EXPECT_EQ(s3_tagQuals_clipped, view.DeletionQVs().Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, view.LabelQVs().Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, view.AltLabelQVs().Fastq());
    EXPECT_EQ(s3_frames_clipped,   view.IPD().Data());
}

TEST(BamRecordTest, ClipEncodedFrames)
{
    const Position qStart  = 500;
    const Position qEnd    = 510;
    const string seq       = "AACCGTTAGC";
    const string quals     = "?]?]?]?]?*";
    const string tagBases  = "AACCGTTAGC";
    const string tagQuals  = "?]?]?]?]?*";
    const f_data frames    = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const string pulseCall   = "ttAaAtaCCGggatTTAcatGCt";
    const string pulseBases  = pulseCall;
    const string pulseQuals  = "==?=]==?]?====]?]===?*=";
    const f_data pulseFrames = { 0,0,10,0,10,0,0,20,20,30,0,0,0,0,40,40,10,0,0,0,30,20,0 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const string seq_clipped      = "CCGTTAG";
    const string quals_clipped    = "?]?]?]?";
    const string tagBases_clipped = "CCGTTAG";
    const string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const string pulseCall_clipped = "CCGggatTTAcatG";
    const string pulseQuals_clipped = "?]?====]?]===?";
    const f_data pulseFrames_clipped = { 20,20,30,0,0,0,0,40,40,10,0,0,0,30 };

    const string seq_rev       = "GCTAACGGTT";
    const string pulseCall_rev = "aGCatgTAAatccCGGtaTtTaa";
    const string quals_rev     = "*?]?]?]?]?";
    const string tagQuals_rev  = quals_rev;
    const f_data frames_rev    = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const string seq_rev_clipped   = "CTAACGG";
    const string quals_rev_clipped = "?]?]?]?";
    const string tagBases_rev_clipped = seq_rev_clipped;
    const string tagQuals_rev_clipped = quals_rev_clipped;
    const f_data frames_rev_clipped = { 30, 10, 40, 40, 30, 20, 20 };

    const string pulseCall_rev_clipped = "CatgTAAatccCGG";
    const string pulseQuals_rev_clipped    = "?===]?]====?]?";
    const f_data pulseFrames_rev_clipped = { 30,0,0,0,10,40,40,0,0,0,0,30,20,20 };

    const string s1_cigar = "10=";
    const string s2_cigar = "5=3D5=";
    const string s3_cigar = "4=1D2I2D4=";

    const string s1_cigar_clipped = "7=";
    const string s2_cigar_clipped = "3=3D4=";
    const string s3_cigar_clipped = "2=1D2I2D3=";

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  pulseCall, pulseBases, pulseQuals, pulseFrames, FrameCodec::V1);

    BamRecord s0 = prototype; // unmapped record
    BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
    BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
    BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
    BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
    BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
    BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);

    s0.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s1.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    {   // s0

        EXPECT_FALSE(s0.IsMapped());
        EXPECT_EQ(clipStart, s0.QueryStart());
        EXPECT_EQ(clipEnd,   s0.QueryEnd());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.AlignedStart());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.AlignedEnd());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.ReferenceStart());
        EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.ReferenceEnd());

        const BamRecordView view
        {
            s0,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,       view.Sequence());
        EXPECT_EQ(quals_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_clipped, view.PulseCalls());
    }

    {   // s1 - FORWARD

        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(clipStart, s1.QueryStart());
        EXPECT_EQ(clipEnd,   s1.QueryEnd());
        EXPECT_EQ(clipStart, s1.AlignedStart());   // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s1.AlignedEnd());     // alignStart + seqLength
        EXPECT_EQ(102, s1.ReferenceStart());       // 100 + startOffset
        EXPECT_EQ(109, s1.ReferenceEnd());         // RefStart + 7=

        EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

        const BamRecordView view
        {
            s1,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,       view.Sequence());
        EXPECT_EQ(quals_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_clipped, view.PulseCalls());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s1_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s1_rev.QueryEnd());
        EXPECT_EQ(clipStart, s1_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s1_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(102, s1_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(109, s1_rev.ReferenceEnd());          // RefStart + 7=

        EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_rev_clipped,       view.Sequence());
        EXPECT_EQ(quals_rev_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_rev_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_rev_clipped, view.PulseCalls());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(clipStart, s2.QueryStart());
        EXPECT_EQ(clipEnd,   s2.QueryEnd());
        EXPECT_EQ(clipStart, s2.AlignedStart());   // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s2.AlignedEnd());     // alignStart + seqLength
        EXPECT_EQ(102, s2.ReferenceStart());       // 100 + startOffset
        EXPECT_EQ(112, s2.ReferenceEnd());         // RefStart + 7= + 3D

        EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

        const BamRecordView view
        {
            s2,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,      view.Sequence());
        EXPECT_EQ(quals_clipped,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped, view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped, view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s2_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
        EXPECT_EQ(clipStart, s2_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s2_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(102, s2_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(112, s2_rev.ReferenceEnd());          // RefStart + 7= + 3D

        EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_rev_clipped,       view.Sequence());
        EXPECT_EQ(quals_rev_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_rev_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_rev_clipped, view.PulseCalls());
    }

    {   // s3 - FORWARD

        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(clipStart, s3.QueryStart());
        EXPECT_EQ(clipEnd,   s3.QueryEnd());
        EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
        EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5= + 3D

        EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

        const BamRecordView view
        {
            s3,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_clipped,       view.Sequence());
        EXPECT_EQ(quals_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_clipped, view.PulseCalls());
    }

    {   // s3 - REVERSE

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s3_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
        EXPECT_EQ(clipStart, s3_rev.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s3_rev.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(102, s3_rev.ReferenceStart());         // 100 + startOffset
        EXPECT_EQ(110, s3_rev.ReferenceEnd());           // RefStart + 5= + 3D

        EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

        const BamRecordView view
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq_rev_clipped,       view.Sequence());
        EXPECT_EQ(quals_rev_clipped,     view.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev_clipped,  view.DeletionTags());
        EXPECT_EQ(tagQuals_rev_clipped,  view.DeletionQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.LabelQVs().Fastq());
        EXPECT_EQ(pulseQuals_rev_clipped,  view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev_clipped,    view.IPD().Data());
        EXPECT_EQ(pulseCall_rev_clipped, view.PulseCalls());
    }

}

// clang-format on

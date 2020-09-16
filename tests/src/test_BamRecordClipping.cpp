// Author: Derek Barnett

#include <pbbam/BamRecord.h>

#include <cstdint>

#include <chrono>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamRecordView.h>
#include <pbbam/BamTagCodec.h>
#include <pbbam/EntireFileQuery.h>

#include "PbbamTestData.h"

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;

using f_data = std::vector<uint16_t>;

namespace BamRecordClippingTests {

static
ReadGroupInfo MakeReadGroup(const Data::FrameCodec codec,
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
                     const std::string& seq,
                     const std::string& quals,
                     const std::string& tagBases,
                     const std::string& tagQuals,
                     const f_data& frames,
                     const std::string& pulseCall = "",
                     const std::string& pulseBases = "",
                     const std::string& pulseQuals = "",
                     const f_data& pulseFrames = f_data(),
                     const Data::FrameCodec codec = Data::FrameCodec::RAW)
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
BamRecord MakeCCSRecord(const std::string& seq,
                        const std::string& quals,
                        const std::string& tagBases,
                        const std::string& tagQuals,
                        const f_data& frames,
                        const std::string& pulseCall = "",
                        const std::string& pulseBases = "",
                        const std::string& pulseQuals = "",
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

TEST(BAM_BamRecordClipping, correctly_performs_clip_to_query_simple)
{
    const Position qStart  = 500;
    const Position qEnd    = 510;
    const std::string seq       = "AACCGTTAGC";
    const std::string quals     = "?]?]?]?]?*";
    const std::string tagBases  = "AACCGTTAGC";
    const std::string tagQuals  = "?]?]?]?]?*";
    const f_data frames    = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const std::string pulseCall   = "ttAaAtaCCGggatTTAcatGCt";
    const std::string& pulseBases  = pulseCall;
    const std::string pulseQuals  = "==?=]==?]?====]?]===?*=";
    const f_data pulseFrames = { 0,0,10,0,10,0,0,20,20,30,0,0,0,0,40,40,10,0,0,0,30,20,0 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const std::string seq_clipped      = "CCGTTAG";
    const std::string quals_clipped    = "?]?]?]?";
    const std::string tagBases_clipped = "CCGTTAG";
    const std::string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const std::string pulseCall_clipped = "CCGggatTTAcatG";
    const std::string pulseQuals_clipped = "?]?====]?]===?";
    const f_data pulseFrames_clipped = { 20,20,30,0,0,0,0,40,40,10,0,0,0,30 };

    const std::string seq_rev       = "GCTAACGGTT";
    const std::string pulseCall_rev = "aGCatgTAAatccCGGtaTtTaa";
    const std::string quals_rev     = "*?]?]?]?]?";
    const f_data frames_rev    = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const std::string seq_rev_clipped   = "CTAACGG";
    const std::string quals_rev_clipped = "?]?]?]?";
    const std::string& tagBases_rev_clipped = seq_rev_clipped;
    const std::string& tagQuals_rev_clipped = quals_rev_clipped;
    const f_data frames_rev_clipped = { 30, 10, 40, 40, 30, 20, 20 };

    const std::string pulseCall_rev_clipped = "CatgTAAatccCGG";
    const std::string pulseQuals_rev_clipped    = "?===]?]====?]?";
    const f_data pulseFrames_rev_clipped = { 30,0,0,0,10,40,40,0,0,0,0,30,20,20 };

    const std::string s1_cigar = "10=";
    const std::string s2_cigar = "5=3D5=";
    const std::string s3_cigar = "4=1D2I2D4=";

    const std::string s1_cigar_clipped = "7=";
    const std::string s2_cigar_clipped = "3=3D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D3=";

    const std::string s1_rev_cigar_clipped = "7=";
    const std::string s2_rev_cigar_clipped = "4=3D3=";
    const std::string s3_rev_cigar_clipped = "3=1D2I2D2=";

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  pulseCall, pulseBases, pulseQuals, pulseFrames);

    {
        SCOPED_TRACE("s0");

        BamRecord s0 = prototype; // unmapped record
        s0.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s1 - FORWARD");

        BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
        s1.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s1 - REVERSE");

        BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(500,  s1_rev.QueryStart());
        EXPECT_EQ(510,  s1_rev.QueryEnd());
        EXPECT_EQ(500,  s1_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(510,  s1_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos, s1_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(110,  s1_rev.ReferenceEnd());          // RefStart + 7=
        EXPECT_EQ(s1_cigar, s1_rev.CigarData().ToStdString());

        s1_rev.Clip(ClipType::CLIP_TO_QUERY, 502, 509);
        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(502, s1_rev.QueryStart());
        EXPECT_EQ(509, s1_rev.QueryEnd());
        EXPECT_EQ(502, s1_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(509, s1_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(102, s1_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(109, s1_rev.ReferenceEnd());          // RefStart + 7=
        EXPECT_EQ(s1_rev_cigar_clipped, s1_rev.CigarData().ToStdString());

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
    {
        SCOPED_TRACE("s2 - FORWARD");

        BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
        s2.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s2 - REVERSE");

        BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
        s2_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s2_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
        EXPECT_EQ(clipStart, s2_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s2_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(102, s2_rev.ReferenceStart());        // 100 + startOffset
        EXPECT_EQ(112, s2_rev.ReferenceEnd());          // RefStart + 7= + 3D

        EXPECT_EQ(s2_rev_cigar_clipped, s2_rev.CigarData().ToStdString());

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
    {
        SCOPED_TRACE("s3 - FORWARD");

        BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
        s3.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s3 - REVERSE");

        BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);
        s3_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s3_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
        EXPECT_EQ(clipStart, s3_rev.AlignedStart());     // queryStart (no soft clips)
        EXPECT_EQ(clipEnd,   s3_rev.AlignedEnd());       // alignStart + seqLength
        EXPECT_EQ(102, s3_rev.ReferenceStart());         // 100 + startOffset
        EXPECT_EQ(110, s3_rev.ReferenceEnd());           // RefStart + 5= + 3D

        EXPECT_EQ(s3_rev_cigar_clipped, s3_rev.CigarData().ToStdString());

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

TEST(BAM_BamRecordClipping, correctly_performs_clip_to_query_with_soft_clips)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const std::string seq      = "TTAACCGTTAGCAAA";
    const std::string seq_rev  = "TTTGCTAACGGTTAA";
    const std::string quals    = "--?]?]?]?]?*+++";
    const std::string tagBases = "TTAACCGTTAGCAAA";
    const std::string tagQuals = "--?]?]?]?]?*+++";
    const std::string tagQuals_rev = "+++*?]?]?]?]?--";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };
    const f_data frames_rev = { 10, 10, 10, 20, 30, 10, 40, 40, 30, 20, 20, 10, 10, 40, 40 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const std::string s1_cigar = "2S10=3S";
    const std::string s1_cigar_clipped = "7=";
    const std::string s1_seq_clipped      = "AACCGTT";
    const std::string s1_quals_clipped    = "?]?]?]?";
    const std::string& s1_tagBases_clipped = s1_seq_clipped;
    const std::string& s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const std::string s1_cigar_rev_clipped = "6=1S";
    const std::string s1_seq_rev_clipped   = "AACGGTT";
    const std::string s1_quals_rev_clipped = "?]?]?]?";
    const std::string& s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const std::string& s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 40, 40, 30, 20, 20, 10, 10 };

    const std::string s2_cigar = "2S5=3D5=3S";
    const std::string s2_cigar_clipped = "5=3D2=";
    const std::string s2_seq_clipped      = "AACCGTT";
    const std::string s2_quals_clipped    = "?]?]?]?";
    const std::string& s2_tagBases_clipped = s2_seq_clipped;
    const std::string& s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const std::string s2_cigar_rev_clipped = "1=3D5=1S";
    const std::string s2_seq_rev_clipped   = "AACGGTT";
    const std::string s2_quals_rev_clipped = "?]?]?]?";
    const std::string& s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const std::string& s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 40, 40, 30, 20, 20, 10, 10 };

    const std::string s3_cigar = "2S4=1D2I2D4=3S";
    const std::string s3_cigar_clipped = "4=1D2I2D1=";
    const std::string s3_seq_clipped      = "AACCGTT";
    const std::string s3_quals_clipped    = "?]?]?]?";
    const std::string& s3_tagBases_clipped = s3_seq_clipped;
    const std::string& s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const std::string s3_cigar_rev_clipped = "1D2I2D4=1S";
    const std::string s3_seq_rev_clipped   = "AACGGTT";
    const std::string s3_quals_rev_clipped = "?]?]?]?";
    const std::string& s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const std::string& s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
    const f_data s3_frames_rev_clipped = { 40, 40, 30, 20, 20, 10, 10 };

    const BamRecord prototype = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                                  seq, tagBases, tagQuals, frames);

    {
        SCOPED_TRACE("s1 - FORWARD");

        BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(100, s1.ReferenceStart());
        EXPECT_EQ(110, s1.ReferenceEnd()); // 10=

        s1.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s1 - REVERSE");

        BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(100, s1_rev.ReferenceStart());
        EXPECT_EQ(110, s1_rev.ReferenceEnd()); // 10=

        s1_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s1_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s1_rev.QueryEnd());
        EXPECT_EQ(503, s1_rev.AlignedStart());    // queryStart (no soft clips)
        EXPECT_EQ(509,   s1_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s1_rev.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 6,  s1_rev.ReferenceEnd());    // RefStart + 7=

        EXPECT_EQ(s1_cigar_rev_clipped, s1_rev.CigarData().ToStdString());

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
    {
        SCOPED_TRACE("s2 - FORWARD");

        BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(100, s2.ReferenceStart());
        EXPECT_EQ(113, s2.ReferenceEnd());   // 5= + 3D + 5=

        s2.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s2 - REVERSE");

        BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(100, s2_rev.ReferenceStart());
        EXPECT_EQ(113, s2_rev.ReferenceEnd());   // 5= + 3D + 5=

        s2_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s2_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
        EXPECT_EQ(503, s2_rev.AlignedStart());    // queryStart (no soft clips left)
        EXPECT_EQ(509,   s2_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s2_rev.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 9, s2_rev.ReferenceEnd());    // RefStart + 5=3D2=

        EXPECT_EQ(s2_cigar_rev_clipped, s2_rev.CigarData().ToStdString());

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
    {
        SCOPED_TRACE("s3 - FORWARD");

        BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(100, s3.ReferenceStart());
        EXPECT_EQ(111, s3.ReferenceEnd());   // 4= + 1D + 2D + 4=

        s3.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

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
    {
        SCOPED_TRACE("s3 - REVERSE");

        BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);
        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(100, s3_rev.ReferenceStart());
        EXPECT_EQ(111, s3_rev.ReferenceEnd());   // 4= + 1D + 2D + 4=

        s3_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(clipStart, s3_rev.QueryStart());
        EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
        EXPECT_EQ(503, s3_rev.AlignedStart());    // queryStart + 1S
        EXPECT_EQ(509,   s3_rev.AlignedEnd());      // alignStart + seqLength
        EXPECT_EQ(tPos,      s3_rev.ReferenceStart());  // 100 + startOffset
        EXPECT_EQ(tPos + 7,  s3_rev.ReferenceEnd());    // RefStart + 4=1D2D1=

        EXPECT_EQ(s3_cigar_rev_clipped, s3_rev.CigarData().ToStdString());

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

TEST(BAM_BamRecordClipping, correctly_performs_clip_to_reference_simple)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const std::string tagQuals_rev = "*?]?]?]?]?";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const std::string s1_cigar = "10=";
    const std::string s1_cigar_clipped = "5=";
    const std::string s1_seq_clipped      = "CCGTT";
    const std::string s1_quals_clipped    = "?]?]?";
    const std::string& s1_tagBases_clipped = s1_seq_clipped;
    const std::string& s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 20, 20, 30, 40, 40 };
    const std::string s1_seq_rev_clipped   = "TAACG";
    const std::string s1_quals_rev_clipped = "]?]?]";
    const std::string& s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const std::string& s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 10, 40, 40, 30, 20 };

    const std::string s2_cigar = "5=3D5=";
    const std::string s2_cigar_clipped = "3=2D";
    const std::string s2_seq_clipped      = "CCG";
    const std::string s2_quals_clipped    = "?]?";
    const std::string& s2_tagBases_clipped = s2_seq_clipped;
    const std::string& s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 20, 20, 30 };
    const std::string s2_seq_rev_clipped   = "TAA";
    const std::string s2_quals_rev_clipped = "]?]";
    const std::string& s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const std::string& s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 10, 40, 40 };

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D";
    const std::string s3_seq_clipped      = "CCGT";
    const std::string s3_quals_clipped    = "?]?]";
    const std::string& s3_tagBases_clipped = s3_seq_clipped;
    const std::string& s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };
    const std::string s3_seq_rev_clipped   = "TAAC";
    const std::string s3_quals_rev_clipped = "]?]?";
    const std::string& s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const std::string& s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
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

TEST(BAM_BamRecordClipping, correctly_performs_clip_to_reference_with_soft_clips)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const std::string seq      = "TTAACCGTTAGCAAA";
    const std::string quals    = "--?]?]?]?]?*+++";
    const std::string tagBases = "TTAACCGTTAGCAAA";
    const std::string tagQuals = "--?]?]?]?]?*+++";
    const std::string tagQuals_rev = "+++*?]?]?]?]?--";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const std::string seq_rev      = "TTTGCTAACGGTTAA";
    const std::string quals_rev    = "+++*?]?]?]?]?--";
    const f_data frames_rev   = { 10, 10, 10, 20, 30, 10, 40, 40, 30, 20, 20, 10, 10, 40, 40 };

    const std::string s1_cigar = "2S10=3S";
    const std::string s1_cigar_clipped = "5=";
    const std::string s1_seq_clipped      = "CCGTT";
    const std::string s1_quals_clipped    = "?]?]?";
    const std::string& s1_tagBases_clipped = s1_seq_clipped;
    const std::string& s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 20, 20, 30, 40, 40 };
    const std::string s1_seq_rev_clipped   = "CTAAC";
    const std::string s1_quals_rev_clipped = "?]?]?";
    const std::string& s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const std::string& s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 30, 10, 40, 40, 30 };

    const std::string s2_cigar = "2S5=3D5=3S";
    const std::string s2_cigar_clipped = "3=2D";
    const std::string s2_seq_clipped      = "CCG";
    const std::string s2_quals_clipped    = "?]?";
    const std::string& s2_tagBases_clipped = s2_seq_clipped;
    const std::string& s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 20, 20, 30 };
    const std::string s2_seq_rev_clipped   = "CTA";
    const std::string s2_quals_rev_clipped = "?]?";
    const std::string& s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const std::string& s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 30, 10, 40 };

    const std::string s3_cigar = "2S4=1D2I2D4=3S";
    const std::string s3_cigar_clipped = "2=1D2I2D";
    const std::string s3_seq_clipped      = "CCGT";
    const std::string s3_quals_clipped    = "?]?]";
    const std::string& s3_tagBases_clipped = s3_seq_clipped;
    const std::string& s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };
    const std::string s3_seq_rev_clipped   = "CTAA";
    const std::string s3_quals_rev_clipped = "?]?]";
    const std::string& s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const std::string& s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
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

TEST(BAM_BamRecordClipping, can_create_new_record_clipped_to_query)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const std::string seq_clipped      = "CCGTTAG";
    const std::string quals_clipped    = "?]?]?]?";
    const std::string tagBases_clipped = "CCGTTAG";
    const std::string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D3=";

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

TEST(BAM_BamRecordClipping, can_create_new_record_clipped_to_reference)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D";
    const std::string s3_seq_clipped      = "CCGT";
    const std::string s3_quals_clipped    = "?]?]";
    const std::string& s3_tagBases_clipped = s3_seq_clipped;
    const std::string& s3_tagQuals_clipped = s3_quals_clipped;
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

TEST(BAM_BamRecordClipping, can_create_new_record_clipped_to_query_static_method)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const std::string seq_clipped      = "CCGTTAG";
    const std::string quals_clipped    = "?]?]?]?";
    const std::string tagBases_clipped = "CCGTTAG";
    const std::string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D3=";

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

TEST(BAM_BamRecordClipping, can_create_new_record_clipped_to_reference_static_method)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D";
    const std::string s3_seq_clipped      = "CCGT";
    const std::string s3_quals_clipped    = "?]?]";
    const std::string& s3_tagBases_clipped = s3_seq_clipped;
    const std::string& s3_tagQuals_clipped = s3_quals_clipped;
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

TEST(BAM_BamRecordClipping, correctly_clips_cigar)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const std::string seq      = "TTAACCGTTAGCAAA";
    const std::string quals    = "--?]?]?]?]?*+++";
    const std::string tagBases = "TTAACCGTTAGCAAA";
    const std::string tagQuals = "--?]?]?]?]?*+++";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };
    const uint8_t mapQual = 80;
    BamRecord s3 = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                     seq, tagBases, tagQuals, frames);
    BamRecord s3_rev = BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames,
                                         seq, tagBases, tagQuals, frames);

    const std::string s3_cigar = "5H2S4=1D2I2D4=3S7H";
    s3.Map(0, 100, Strand::FORWARD, s3_cigar, mapQual);
    s3_rev.Map(0, 100, Strand::REVERSE, s3_cigar, mapQual);

    const Cigar s3_cigar_raw     = s3.CigarData();
    const Cigar s3_cigar_clipped = s3.CigarData(true);

    EXPECT_EQ(s3_cigar, s3_cigar_raw.ToStdString());
    EXPECT_EQ(std::string("4=1D2I2D4="), s3_cigar_clipped.ToStdString());
}

TEST(BAM_BamRecordClipping, can_make_ccs_record_clipped_to_query)
{
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 2;
    const Position clipEnd   = 9;

    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const std::string seq_clipped      = "CCGTTAG";
    const std::string quals_clipped    = "?]?]?]?";
    const std::string tagBases_clipped = "CCGTTAG";
    const std::string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D3=";

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

TEST(BAM_BamRecordClipping, can_make_ccs_record_clipped_to_reference)
{
    const std::string seq      = "AACCGTTAGC";
    const std::string quals    = "?]?]?]?]?*";
    const std::string tagBases = "AACCGTTAGC";
    const std::string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const std::string s3_cigar = "4=1D2I2D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D";
    const std::string s3_seq_clipped      = "CCGT";
    const std::string s3_quals_clipped    = "?]?]";
    const std::string& s3_tagBases_clipped = s3_seq_clipped;
    const std::string& s3_tagQuals_clipped = s3_quals_clipped;
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

TEST(BAM_BamRecordClipping, correctly_clips_encoded_frames)
{
    const Position qStart  = 500;
    const Position qEnd    = 510;
    const std::string seq       = "AACCGTTAGC";
    const std::string quals     = "?]?]?]?]?*";
    const std::string tagBases  = "AACCGTTAGC";
    const std::string tagQuals  = "?]?]?]?]?*";
    const f_data frames    = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const std::string pulseCall   = "ttAaAtaCCGggatTTAcatGCt";
    const std::string& pulseBases  = pulseCall;
    const std::string pulseQuals  = "==?=]==?]?====]?]===?*=";
    const f_data pulseFrames = { 0,0,10,0,10,0,0,20,20,30,0,0,0,0,40,40,10,0,0,0,30,20,0 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const std::string seq_clipped      = "CCGTTAG";
    const std::string quals_clipped    = "?]?]?]?";
    const std::string tagBases_clipped = "CCGTTAG";
    const std::string tagQuals_clipped = "?]?]?]?";
    const f_data frames_clipped   = { 20, 20, 30, 40, 40, 10, 30 };

    const std::string pulseCall_clipped = "CCGggatTTAcatG";
    const std::string pulseQuals_clipped = "?]?====]?]===?";
    const f_data pulseFrames_clipped = { 20,20,30,0,0,0,0,40,40,10,0,0,0,30 };

    const std::string seq_rev       = "GCTAACGGTT";
    const std::string pulseCall_rev = "aGCatgTAAatccCGGtaTtTaa";
    const std::string quals_rev     = "*?]?]?]?]?";
    const f_data frames_rev    = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const std::string seq_rev_clipped   = "CTAACGG";
    const std::string quals_rev_clipped = "?]?]?]?";
    const std::string& tagBases_rev_clipped = seq_rev_clipped;
    const std::string& tagQuals_rev_clipped = quals_rev_clipped;
    const f_data frames_rev_clipped = { 30, 10, 40, 40, 30, 20, 20 };

    const std::string pulseCall_rev_clipped = "CatgTAAatccCGG";
    const std::string pulseQuals_rev_clipped    = "?===]?]====?]?";
    const f_data pulseFrames_rev_clipped = { 30,0,0,0,10,40,40,0,0,0,0,30,20,20 };

    const std::string s1_cigar = "10=";
    const std::string s2_cigar = "5=3D5=";
    const std::string s3_cigar = "4=1D2I2D4=";

    const std::string s1_cigar_clipped = "7=";
    const std::string s2_cigar_clipped = "3=3D4=";
    const std::string s3_cigar_clipped = "2=1D2I2D3=";

    const std::string s1_cigar_rev_clipped = "7=";
    const std::string s2_cigar_rev_clipped = "4=3D3=";
    const std::string s3_cigar_rev_clipped = "3=1D2I2D2=";

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

        EXPECT_EQ(s1_cigar_rev_clipped, s1_rev.CigarData().ToStdString());

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

        EXPECT_EQ(s2_cigar_rev_clipped, s2_rev.CigarData().ToStdString());

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

        EXPECT_EQ(s3_cigar_rev_clipped, s3_rev.CigarData().ToStdString());

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

TEST(BAM_BamRecordClipping, can_excise_soft_clips_from_frames_with_deletions)
{
    const std::string expectedName{"m141008_060349_42194_c100704972550000001823137703241586_s1_p0/14/2409_2745"};
    const PacBio::BAM::Strand expectedStrand = PacBio::BAM::Strand::FORWARD;
    const std::string expectedCigar{
        "20S11=1I47=1I2=1I6=1I22=1I2=1I9=1I29=1D6=1I16=1I6=1I7=1I8=2I5=1I5=1I11=1I5=5I2=3I1=1I1=1I1=3I5=2D19=1I14=1I17=28S"};
    const std::string expectedRawSeq{
        "CCCCGGGATTCCTCTAGATGCATCAGGTAAGAAAAGTACGATGCTACAGCTTGTGACTGGTGCGGCACTT"
        "TTGGCTGAGTTTATCCTGTGCCACCTCATGTATTCTGCCCTAGACAGTCGGTCTTGCACGCCATTACTAG"
        "ACCGACAAAATGGAACCGGGGCCCTTAAACCCCGTTCGAAGGCGTAAGCAAGGAAGATAGGGTTTTATGA"
        "AACTCTTCCCAGTCAATAATACCAAAAAAACCCCAACCAAGATCGTGACGGATTGCAGAGCGAATCCTAT"
        "CCGCGCTCGCAATAATTTAGTGTTGATCCAAGCTTGCTGAGGACTAGTAAAGCTTC"};
    const std::string expectedClippedSeq{
        "CATCAGGTAAGAAAAGTACGATGCTACAGCTTGTGACTGGTGCGGCACTTTTGGCTGAGTTTATCCTGTG"
        "CCACCTCATGTATTCTGCCCTAGACAGTCGGTCTTGCACGCCATTACTAGACCGACAAAATGGAACCGGG"
        "GCCCTTAAACCCCGTTCGAAGGCGTAAGCAAGGAAGATAGGGTTTTATGAAACTCTTCCCAGTCAATAAT"
        "ACCAAAAAAACCCCAACCAAGATCGTGACGGATTGCAGAGCGAATCCTATCCGCGCTCGCAATAATTTAG"
        "TGTTGATC"};
    const std::vector<uint8_t> expectedRawIpds{
        17,3,8,3,4,1,14,8,2,1,21,3,1,17,22,13,10,9,89,7,4,5,3,17,8,8,18,58,14,
        25,8,5,9,1,5,0,20,16,15,9,78,19,2,20,23,12,2,5,7,3,5,61,19,12,13,6,65,
        18,105,2,34,94,3,38,69,16,5,76,1,21,5,3,2,0,32,23,26,9,3,4,18,2,2,12,19,
        33,63,11,4,25,3,7,7,3,26,48,28,34,1,2,6,31,17,29,68,5,20,79,6,12,10,3,
        43,72,21,65,8,45,17,14,13,20,7,3,5,8,0,17,11,65,6,7,8,3,6,11,4,1,80,4,
        16,21,12,4,2,8,1,25,22,36,18,34,11,5,4,33,3,12,1,14,8,22,4,8,76,8,5,18,
        32,5,33,47,255,36,9,26,2,6,47,0,35,8,8,0,5,37,40,1,11,8,39,60,8,42,0,3,
        6,11,12,20,24,15,1,10,10,38,25,63,21,28,0,4,17,0,31,23,13,41,23,42,0,7,
        33,7,23,11,50,30,2,44,21,182,44,105,231,33,255,59,189,253,17,13,7,28,40,
        84,8,13,34,70,214,174,103,5,8,1,8,9,8,1,12,7,4,17,7,45,2,2,7,10,7,19,28,
        31,3,18,0,42,0,8,2,9,2,1,11,25,1,35,36,1,7,5,17,12,39,8,31,1,40,41,4,18,
        2,51,14,1,16,255,2,5,83,2,6,2,1,6,9,10,3,31,19,35,6,16,21,12,28,4,10,10,
        12,1,105,17,2,11};
    const std::vector<uint8_t> expectedClippedIpds{
        4,5,3,17,8,8,18,58,14,25,8,5,9,1,5,0,20,16,15,9,78,19,2,20,23,12,2,5,7,
        3,5,61,19,12,13,6,65,18,105,2,34,94,3,38,69,16,5,76,1,21,5,3,2,0,32,23,
        26,9,3,4,18,2,2,12,19,33,63,11,4,25,3,7,7,3,26,48,28,34,1,2,6,31,17,29,
        68,5,20,79,6,12,10,3,43,72,21,65,8,45,17,14,13,20,7,3,5,8,0,17,11,65,6,
        7,8,3,6,11,4,1,80,4,16,21,12,4,2,8,1,25,22,36,18,34,11,5,4,33,3,12,1,14,
        8,22,4,8,76,8,5,18,32,5,33,47,255,36,9,26,2,6,47,0,35,8,8,0,5,37,40,1,
        11,8,39,60,8,42,0,3,6,11,12,20,24,15,1,10,10,38,25,63,21,28,0,4,17,0,31,
        23,13,41,23,42,0,7,33,7,23,11,50,30,2,44,21,182,44,105,231,33,255,59,
        189,253,17,13,7,28,40,84,8,13,34,70,214,174,103,5,8,1,8,9,8,1,12,7,4,17,
        7,45,2,2,7,10,7,19,28,31,3,18,0,42,0,8,2,9,2,1,11,25,1,35,36,1,7,5,17,
        12,39,8,31,1,40,41,4,18,2,51,14,1,16,255};

    const std::string fn{PbbamTestsConfig::Data_Dir + "/softclip_deletions.bam"};
    BamRecord record;
    BamReader reader{fn};
    ASSERT_TRUE(reader.GetNext(record));

    EXPECT_EQ(expectedName, record.FullName());
    EXPECT_EQ(expectedStrand, record.AlignedStrand());
    EXPECT_EQ(expectedCigar, record.CigarData().ToStdString());

    const auto rawSeq = record.Sequence(PacBio::BAM::Orientation::GENOMIC);
    const auto clippedSeq = record.Sequence(PacBio::BAM::Orientation::GENOMIC, false, true);
    EXPECT_EQ(expectedRawSeq, rawSeq);
    EXPECT_EQ(expectedClippedSeq, clippedSeq);

    ASSERT_TRUE(record.HasIPD());
    const auto rawIpds = record.IPD(PacBio::BAM::Orientation::GENOMIC).Encode();
    const auto clippedIpds = record.IPD(PacBio::BAM::Orientation::GENOMIC, false, true).Encode();
    EXPECT_EQ(expectedRawIpds, rawIpds);
    EXPECT_EQ(expectedClippedIpds, clippedIpds);
}

TEST(BAM_BamRecordClipping, can_clip_to_query_stranded)
{
    using namespace PacBio::BAM;

    const std::string bamFile{PbbamTestsConfig::Data_Dir + "/clip_to_query.bam"};

    bool first = true;
    EntireFileQuery query{bamFile};
    for (auto& i : query)
    {
        Strand expectedStrand;
        std::string scope;
        if (first) {
            expectedStrand = Strand::FORWARD;
            scope = "First record (FORWARD strand)";
        } else {
            expectedStrand = Strand::REVERSE;
            scope = "Second record (REVERSE strand)";
        }

        SCOPED_TRACE(scope);

        // initial
        EXPECT_EQ(2, i.ReferenceStart());
        EXPECT_EQ(7, i.ReferenceEnd());
        EXPECT_EQ(0, i.QueryStart());
        EXPECT_EQ(8, i.QueryEnd());
        EXPECT_EQ(expectedStrand, i.AlignedStrand());
        EXPECT_EQ("1S4=1I1=1S", i.CigarData().ToStdString());

        // first clip
        i.Clip(ClipType::CLIP_TO_REFERENCE, 3, 6);
        EXPECT_EQ(3, i.ReferenceStart());
        EXPECT_EQ(6, i.ReferenceEnd());
        EXPECT_EQ(2, i.QueryStart());
        EXPECT_EQ(6, i.QueryEnd());
        EXPECT_EQ(expectedStrand, i.AlignedStrand());
        EXPECT_EQ("3=1I", i.CigarData().ToStdString());

        // second clip
        Position qS;
        Position qE;
        if (first) {
            qS = i.QueryStart();
            qE = i.QueryEnd() - 1;
        } else {
            qS = i.QueryStart() + 1;
            qE = i.QueryEnd();
        }
        i.Clip(ClipType::CLIP_TO_QUERY, qS, qE);
        EXPECT_EQ(3, i.ReferenceStart());
        EXPECT_EQ(6, i.ReferenceEnd());
        EXPECT_EQ(qS, i.QueryStart());
        EXPECT_EQ(qE, i.QueryEnd());
        EXPECT_EQ(expectedStrand, i.AlignedStrand());
        EXPECT_EQ("3=", i.CigarData().ToStdString());

        first = false;
    }
}

TEST(BAM_BamRecordClipping, clipping_flanking_inserts_is_ignored_on_clip_to_query)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const std::string seq      = "TTAACCGTTAGCAAA";
    const std::string quals    = "--?]?]?]?]?*+++";
    const std::string tagBases = "TTAACCGTTAGCAAA";
    const std::string tagQuals = "--?]?]?]?]?*+++";
    const f_data frames = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };

    const BamRecord prototype =
        BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals,
                                           frames, seq, tagBases, tagQuals, frames);

    {   // aligned forward

        const int32_t  tId     = 0;
        const Position tPos    = 100;
        const uint8_t  mapQual = 80;
        const Cigar cigar{"4I5=6I"};

        BamRecord s = prototype.Mapped(tId, tPos, Strand::FORWARD, cigar, mapQual);
        EXPECT_TRUE(s.IsMapped());
        EXPECT_EQ(100, s.ReferenceStart());
        EXPECT_EQ(105, s.ReferenceEnd());

        const size_t clipStart = 502;
        const size_t clipEnd = 512;
        const bool exciseFlankingInserts = true;

        s.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd, exciseFlankingInserts);

        EXPECT_TRUE(s.IsMapped());
        EXPECT_EQ(Strand::FORWARD, s.AlignedStrand());
        EXPECT_EQ("2I5=3I", s.CigarData().ToStdString());

        EXPECT_EQ(clipStart, s.QueryStart());
        EXPECT_EQ(clipEnd,   s.QueryEnd());
        EXPECT_EQ(clipStart, s.AlignedStart());
        EXPECT_EQ(clipEnd,   s.AlignedEnd());
        EXPECT_EQ(100,       s.ReferenceStart());
        EXPECT_EQ(105,       s.ReferenceEnd());
    }
    {   // aligned reverse

        const int32_t  tId     = 0;
        const Position tPos    = 100;
        const uint8_t  mapQual = 80;
        const Cigar cigar{"4I5=6I"};

        BamRecord s = prototype.Mapped(tId, tPos, Strand::REVERSE, cigar, mapQual);
        EXPECT_TRUE(s.IsMapped());
        EXPECT_EQ(100, s.ReferenceStart());
        EXPECT_EQ(105, s.ReferenceEnd());

        const size_t clipStart = 502;
        const size_t clipEnd = 512;
        const bool exciseFlankingInserts = true;

        s.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd, exciseFlankingInserts);

        EXPECT_TRUE(s.IsMapped());
        EXPECT_EQ(Strand::REVERSE, s.AlignedStrand());
        EXPECT_EQ("1I5=4I", s.CigarData().ToStdString());

        EXPECT_EQ(clipStart, s.QueryStart());
        EXPECT_EQ(clipEnd,   s.QueryEnd());
        EXPECT_EQ(clipStart, s.AlignedStart());
        EXPECT_EQ(clipEnd,   s.AlignedEnd());
        EXPECT_EQ(100,       s.ReferenceStart());
        EXPECT_EQ(105,       s.ReferenceEnd());
    }
}

TEST(BAM_BamRecordClipping, can_excise_flanking_insertsion_when_clipping_to_reference_forward)
{
    const Position qStart = 500;
    const Position qEnd   = 526;
    const std::string seq      = "TTAACCGTTAGCAAATTAACCGTTAG";
    const std::string quals    = "--?]?]?]?]?*+++--?]?]?]?]?";
    const std::string tagBases = "TTAACCGTTAGCAAATTAACCGTTAG";
    const std::string tagQuals = "--?]?]?]?]?*+++--?]?]?]?]?";
    const f_data frames = {
        40, 40, 10, 10, 20, 20, 30, 40, 40, 10,
        30, 20, 10, 10, 10, 40, 40, 10, 10, 20,
        20, 30, 40, 40, 10, 30 };

    const BamRecord prototype =
        BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals,
                                           frames, seq, tagBases, tagQuals, frames);

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Cigar cigar{"3=6I10=6I1="};

    const size_t clipStart = 103;
    const size_t clipEnd = 113;

    // ----------------
    // keep inserts

    bool exciseFlankingInserts = false;

    BamRecord withInserts = prototype.Mapped(tId, tPos, Strand::FORWARD, cigar, mapQual);
    EXPECT_TRUE(withInserts.IsMapped());
    EXPECT_EQ(100, withInserts.ReferenceStart());
    EXPECT_EQ(114, withInserts.ReferenceEnd());

    withInserts.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd, exciseFlankingInserts);

    EXPECT_TRUE(withInserts.IsMapped());
    EXPECT_EQ(Strand::FORWARD, withInserts.AlignedStrand());
    EXPECT_EQ("6I10=6I", withInserts.CigarData().ToStdString());

    EXPECT_EQ(503, withInserts.QueryStart());
    EXPECT_EQ(525, withInserts.QueryEnd());
    EXPECT_EQ(503, withInserts.AlignedStart());
    EXPECT_EQ(525, withInserts.AlignedEnd());
    EXPECT_EQ(103, withInserts.ReferenceStart());
    EXPECT_EQ(113, withInserts.ReferenceEnd());

    // -----------------
    // excise inserts

    exciseFlankingInserts = true;

    BamRecord withoutInserts = prototype.Mapped(tId, tPos, Strand::FORWARD, cigar, mapQual);
    EXPECT_TRUE(withoutInserts.IsMapped());
    EXPECT_EQ(100, withoutInserts.ReferenceStart());
    EXPECT_EQ(114, withoutInserts.ReferenceEnd());

    withoutInserts.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd, exciseFlankingInserts);

    EXPECT_TRUE(withoutInserts.IsMapped());
    EXPECT_EQ(Strand::FORWARD, withoutInserts.AlignedStrand());
    EXPECT_EQ("10=", withoutInserts.CigarData().ToStdString());

    EXPECT_EQ(509, withoutInserts.QueryStart());
    EXPECT_EQ(519, withoutInserts.QueryEnd());
    EXPECT_EQ(509, withoutInserts.AlignedStart());
    EXPECT_EQ(519, withoutInserts.AlignedEnd());
    EXPECT_EQ(103, withoutInserts.ReferenceStart());
    EXPECT_EQ(113, withoutInserts.ReferenceEnd());
}

TEST(BAM_BamRecordClipping, can_excise_flanking_insertsion_when_clipping_to_reference_reverse)
{
    const Position qStart = 500;
    const Position qEnd   = 526;
    const std::string seq      = "TTAACCGTTAGCAAATTAACCGTTAG";
    const std::string quals    = "--?]?]?]?]?*+++--?]?]?]?]?";
    const std::string tagBases = "TTAACCGTTAGCAAATTAACCGTTAG";
    const std::string tagQuals = "--?]?]?]?]?*+++--?]?]?]?]?";
    const f_data frames = {
        40, 40, 10, 10, 20, 20, 30, 40, 40, 10,
        30, 20, 10, 10, 10, 40, 40, 10, 10, 20,
        20, 30, 40, 40, 10, 30 };

    const BamRecord prototype =
        BamRecordClippingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals,
                                           frames, seq, tagBases, tagQuals, frames);

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;
    const Cigar cigar{"3=6I10=6I1="};

    const size_t clipStart = 103;
    const size_t clipEnd = 113;

    // ----------------
    // keep inserts

    bool exciseFlankingInserts = false;

    BamRecord withInserts = prototype.Mapped(tId, tPos, Strand::REVERSE, cigar, mapQual);

    EXPECT_TRUE(withInserts.IsMapped());
    EXPECT_EQ(100, withInserts.ReferenceStart());
    EXPECT_EQ(114, withInserts.ReferenceEnd());

    withInserts.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd, exciseFlankingInserts);

    EXPECT_TRUE(withInserts.IsMapped());
    EXPECT_EQ(Strand::REVERSE, withInserts.AlignedStrand());
    EXPECT_EQ("6I10=6I", withInserts.CigarData().ToStdString());

    EXPECT_EQ(501, withInserts.QueryStart());
    EXPECT_EQ(523, withInserts.QueryEnd());
    EXPECT_EQ(501, withInserts.AlignedStart());
    EXPECT_EQ(523, withInserts.AlignedEnd());
    EXPECT_EQ(103, withInserts.ReferenceStart());
    EXPECT_EQ(113, withInserts.ReferenceEnd());

    // -----------------
    // excise inserts

    exciseFlankingInserts = true;

    BamRecord withoutInserts = prototype.Mapped(tId, tPos, Strand::REVERSE, cigar, mapQual);
    EXPECT_TRUE(withoutInserts.IsMapped());
    EXPECT_EQ(100, withoutInserts.ReferenceStart());
    EXPECT_EQ(114, withoutInserts.ReferenceEnd());

    withoutInserts.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd, exciseFlankingInserts);

    EXPECT_TRUE(withoutInserts.IsMapped());
    EXPECT_EQ(Strand::REVERSE, withoutInserts.AlignedStrand());
    EXPECT_EQ("10=", withoutInserts.CigarData().ToStdString());

    EXPECT_EQ(507, withoutInserts.QueryStart());
    EXPECT_EQ(517, withoutInserts.QueryEnd());
    EXPECT_EQ(507, withoutInserts.AlignedStart());
    EXPECT_EQ(517, withoutInserts.AlignedEnd());
    EXPECT_EQ(103, withoutInserts.ReferenceStart());
    EXPECT_EQ(113, withoutInserts.ReferenceEnd());

}

// clang-format on

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
#include <pbbam/BamRecord.h>
#include <pbbam/BamTagCodec.h>
#include <chrono>
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

typedef vector<uint16_t> f_data;

namespace tests {

static
BamRecord MakeRecord(const Position qStart,
                     const Position qEnd,
                     const string& seq,
                     const string& quals,
                     const string& tagBases,
                     const string& tagQuals,
                     const f_data& frames,
                     const string& pulseCall = "")
{
    BamRecordImpl impl;
    impl.SetSequenceAndQualities(seq, quals);

    TagCollection tags;
    tags["qs"] = qStart;
    tags["qe"] = qEnd;
    tags["pa"] = frames;
    tags["pm"] = frames;
    tags["ip"] = frames;
    tags["pw"] = frames;
    tags["dt"] = tagBases;
    tags["st"] = tagBases;
    tags["lt"] = tagBases;
    tags["at"] = tagBases;
    tags["dq"] = tagQuals;
    tags["iq"] = tagQuals;
    tags["mq"] = tagQuals;
    tags["sq"] = tagQuals;
    tags["lq"] = tagQuals;
    tags["aq"] = tagQuals;
    tags["pc"] = pulseCall;
    impl.Tags(tags);

    return BamRecord(std::move(impl));
}

} // namespace tests

TEST(BamRecordClippingTest, ClipToQuery_Basic)
{
    const Position qStart  = 500;
    const Position qEnd    = 510;
    const string seq       = "AACCGTTAGC";
    const string pulseCall = "ttAaAtaCCGggatTTAcatGCt";
    const string quals     = "?]?]?]?]?*";
    const string tagBases  = "AACCGTTAGC";
    const string tagQuals  = "?]?]?]?]?*";
    const f_data frames    = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

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

    const string seq_rev       = "GCTAACGGTT";
    const string pulseCall_rev = "aGCatgTAAatccCGGtaTtTaa";
    const string quals_rev     = "*?]?]?]?]?";
    const f_data frames_rev    = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const string seq_rev_clipped   = "TAACGGT";
    const string quals_rev_clipped = "]?]?]?]";
    const string tagBases_rev_clipped = seq_rev_clipped;
    const string tagQuals_rev_clipped = quals_rev_clipped;
    const f_data frames_rev_clipped = { 10, 40, 40, 30, 20, 20, 10 };

    const string s1_cigar = "10M";
    const string s2_cigar = "5M3D5M";
    const string s3_cigar = "4M1D2I2D4M";

    const string s1_cigar_clipped = "7M";
    const string s2_cigar_clipped = "3M3D4M";
    const string s3_cigar_clipped = "2M1D2I2D3M";

    const BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames, pulseCall);
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

    // s0
    EXPECT_FALSE(s0.IsMapped());
    EXPECT_EQ(clipStart, s0.QueryStart());
    EXPECT_EQ(clipEnd,   s0.QueryEnd());
    EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.AlignedStart());
    EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.AlignedEnd());
    EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.ReferenceStart());
    EXPECT_EQ(PacBio::BAM::UnmappedPosition, s0.ReferenceEnd());
    EXPECT_EQ(seq_clipped,      s0.Sequence());
    EXPECT_EQ(quals_clipped,    s0.Qualities().Fastq());
    EXPECT_EQ(tagBases_clipped, s0.DeletionTag());
    EXPECT_EQ(tagQuals_clipped, s0.DeletionQV().Fastq());
    EXPECT_EQ(tagQuals_clipped, s0.LabelQV().Fastq());
    EXPECT_EQ(tagQuals_clipped, s0.AltLabelQV().Fastq());
    EXPECT_EQ(frames_clipped,   s0.IPD().Data());
    EXPECT_EQ(pulseCall,        s0.PulseCall());

    // s1 - FORWARD
    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
    EXPECT_EQ(clipStart, s1.QueryStart());
    EXPECT_EQ(clipEnd,   s1.QueryEnd());
    EXPECT_EQ(clipStart, s1.AlignedStart());   // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s1.AlignedEnd());     // alignStart + seqLength
    EXPECT_EQ(102, s1.ReferenceStart());       // 100 + startOffset
    EXPECT_EQ(109, s1.ReferenceEnd());         // RefStart + 7M

    EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

    EXPECT_EQ(seq_clipped,      s1.Sequence());
    EXPECT_EQ(quals_clipped,    s1.Qualities().Fastq());
    EXPECT_EQ(tagBases_clipped, s1.DeletionTag());
    EXPECT_EQ(tagQuals_clipped, s1.DeletionQV().Fastq());
    EXPECT_EQ(tagQuals_clipped, s1.LabelQV().Fastq());
    EXPECT_EQ(tagQuals_clipped, s1.AltLabelQV().Fastq());
    EXPECT_EQ(frames_clipped,   s1.IPD().Data());
    EXPECT_EQ(pulseCall,        s1.PulseCall());

    // s1 - REVERSE
    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
    EXPECT_EQ(clipStart, s1_rev.QueryStart());
    EXPECT_EQ(clipEnd,   s1_rev.QueryEnd());
    EXPECT_EQ(clipStart, s1_rev.AlignedStart());    // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s1_rev.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(102, s1_rev.ReferenceStart());        // 100 + startOffset
    EXPECT_EQ(109, s1_rev.ReferenceEnd());          // RefStart + 7M

    EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

    EXPECT_EQ(seq_rev_clipped,      s1_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_rev_clipped,    s1_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_rev_clipped, s1_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_rev_clipped, s1_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_rev_clipped, s1_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_rev_clipped, s1_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_rev_clipped,   s1_rev.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(pulseCall_rev,        s1_rev.PulseCall(Orientation::GENOMIC));

    // s2 - FORWARD
    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
    EXPECT_EQ(clipStart, s2.QueryStart());
    EXPECT_EQ(clipEnd,   s2.QueryEnd());
    EXPECT_EQ(clipStart, s2.AlignedStart());   // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s2.AlignedEnd());     // alignStart + seqLength
    EXPECT_EQ(102, s2.ReferenceStart());       // 100 + startOffset
    EXPECT_EQ(112, s2.ReferenceEnd());         // RefStart + 7M + 3D

    EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

    EXPECT_EQ(seq_clipped,      s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_clipped,    s2.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_clipped, s2.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_clipped, s2.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s2.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s2.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_clipped,   s2.IPD(Orientation::GENOMIC).Data());

    // s2 - REVERSE
    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
    EXPECT_EQ(clipStart, s2_rev.QueryStart());
    EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
    EXPECT_EQ(clipStart, s2_rev.AlignedStart());    // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s2_rev.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(102, s2_rev.ReferenceStart());        // 100 + startOffset
    EXPECT_EQ(112, s2_rev.ReferenceEnd());          // RefStart + 7M + 3D

    EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

    EXPECT_EQ(seq_rev_clipped,      s2_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_rev_clipped,    s2_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_rev_clipped, s2_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_rev_clipped, s2_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_rev_clipped, s2_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_rev_clipped, s2_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_rev_clipped,   s2_rev.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(pulseCall_rev,        s2_rev.PulseCall(Orientation::GENOMIC));

    // s3 - FORWARD
    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(clipStart, s3.QueryStart());
    EXPECT_EQ(clipEnd,   s3.QueryEnd());
    EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5M + 3D

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    EXPECT_EQ(seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());

    // s3 - REVERSE
    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
    EXPECT_EQ(clipStart, s3_rev.QueryStart());
    EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
    EXPECT_EQ(clipStart, s3_rev.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s3_rev.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(102, s3_rev.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3_rev.ReferenceEnd());           // RefStart + 5M + 3D

    EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

    EXPECT_EQ(seq_rev_clipped,      s3_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_rev_clipped,    s3_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_rev_clipped, s3_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_rev_clipped, s3_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_rev_clipped, s3_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_rev_clipped, s3_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_rev_clipped,   s3_rev.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(pulseCall_rev,        s3_rev.PulseCall(Orientation::GENOMIC));
}

TEST(BamRecordClippingTest, ClipToQuery_WithSoftClips)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const string seq      = "TTAACCGTTAGCAAA";
    const string quals    = "--?]?]?]?]?*+++";
    const string tagBases = "TTAACCGTTAGCAAA";
    const string tagQuals = "--?]?]?]?]?*+++";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 502;
    const Position clipEnd   = 509;

    const string s1_cigar = "2S10M3S";
    const string s1_cigar_clipped = "7M";
    const string s1_seq_clipped      = "AACCGTT";
    const string s1_quals_clipped    = "?]?]?]?";
    const string s1_tagBases_clipped = s1_seq_clipped;
    const string s1_tagQuals_clipped = s1_quals_clipped;
    const f_data s1_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const string s1_seq_rev_clipped   = "TGCTAAC";
    const string s1_quals_rev_clipped = "+*?]?]?";
    const string s1_tagBases_rev_clipped = s1_seq_rev_clipped;
    const string s1_tagQuals_rev_clipped = s1_quals_rev_clipped;
    const f_data s1_frames_rev_clipped = { 10, 20, 30, 10, 40, 40, 30 };

    const string s2_cigar = "2S5M3D5M3S";
    const string s2_cigar_clipped = "5M3D2M";
    const string s2_seq_clipped      = "AACCGTT";
    const string s2_quals_clipped    = "?]?]?]?";
    const string s2_tagBases_clipped = s2_seq_clipped;
    const string s2_tagQuals_clipped = s2_quals_clipped;
    const f_data s2_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const string s2_seq_rev_clipped   = "TGCTAAC";
    const string s2_quals_rev_clipped = "+*?]?]?";
    const string s2_tagBases_rev_clipped = s2_seq_rev_clipped;
    const string s2_tagQuals_rev_clipped = s2_quals_rev_clipped;
    const f_data s2_frames_rev_clipped = { 10, 20, 30, 10, 40, 40, 30 };

    const string s3_cigar = "2S4M1D2I2D4M3S";
    const string s3_cigar_clipped = "4M1D2I2D1M";
    const string s3_seq_clipped      = "AACCGTT";
    const string s3_quals_clipped    = "?]?]?]?";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 10, 10, 20, 20, 30, 40, 40 };
    const string s3_seq_rev_clipped   = "TGCTAAC";
    const string s3_quals_rev_clipped = "+*?]?]?";
    const string s3_tagBases_rev_clipped = s3_seq_rev_clipped;
    const string s3_tagQuals_rev_clipped = s3_quals_rev_clipped;
    const f_data s3_frames_rev_clipped = { 10, 20, 30, 10, 40, 40, 30 };

    const BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s1 = prototype.Mapped(tId, tPos, Strand::FORWARD, s1_cigar, mapQual);
    BamRecord s2 = prototype.Mapped(tId, tPos, Strand::FORWARD, s2_cigar, mapQual);
    BamRecord s3 = prototype.Mapped(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);
    BamRecord s1_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s1_cigar, mapQual);
    BamRecord s2_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s2_cigar, mapQual);
    BamRecord s3_rev = prototype.Mapped(tId, tPos, Strand::REVERSE, s3_cigar, mapQual);

    // sanity checks before clipping
    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(tPos, s1.ReferenceStart());
    EXPECT_EQ(tPos + 10, s1.ReferenceEnd()); // 10M

    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(tPos, s1_rev.ReferenceStart());
    EXPECT_EQ(tPos + 10, s1_rev.ReferenceEnd()); // 10M

    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(tPos, s2.ReferenceStart());
    EXPECT_EQ(tPos + 13, s2.ReferenceEnd());   // 5M + 3D + 5M

    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(tPos, s2_rev.ReferenceStart());
    EXPECT_EQ(tPos + 13, s2_rev.ReferenceEnd());   // 5M + 3D + 5M

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(tPos, s3.ReferenceStart());
    EXPECT_EQ(tPos + 11, s3.ReferenceEnd());   // 4M + 1D + 2D + 4M

    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(tPos, s3_rev.ReferenceStart());
    EXPECT_EQ(tPos + 11, s3_rev.ReferenceEnd());   // 4M + 1D + 2D + 4M

    s1.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    // s1 - FORWARD
    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
    EXPECT_EQ(clipStart, s1.QueryStart());
    EXPECT_EQ(clipEnd,   s1.QueryEnd());
    EXPECT_EQ(clipStart, s1.AlignedStart());    // queryStart (no soft clips left)
    EXPECT_EQ(clipEnd,   s1.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(tPos,      s1.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 7,  s1.ReferenceEnd());    // RefStart + 7M

    EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

    EXPECT_EQ(s1_seq_clipped,      s1.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s1_quals_clipped,    s1.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagBases_clipped, s1.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s1_tagQuals_clipped, s1.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_clipped, s1.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_clipped, s1.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_frames_clipped,   s1.IPD(Orientation::GENOMIC).Data());

    // s1 - REVERSE
    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
    EXPECT_EQ(clipStart, s1_rev.QueryStart());
    EXPECT_EQ(clipEnd,   s1_rev.QueryEnd());
    EXPECT_EQ(clipStart, s1_rev.AlignedStart());    // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s1_rev.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(tPos,      s1_rev.ReferenceStart());  // 100 + startOffset
    EXPECT_EQ(tPos + 7,  s1_rev.ReferenceEnd());    // RefStart + 7M

    EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

    EXPECT_EQ(s1_seq_rev_clipped,      s1_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s1_quals_rev_clipped,    s1_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagBases_rev_clipped, s1_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_frames_rev_clipped,   s1_rev.IPD(Orientation::GENOMIC).Data());

    // s2 - FORWARD
    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
    EXPECT_EQ(clipStart, s2.QueryStart());
    EXPECT_EQ(clipEnd,   s2.QueryEnd());
    EXPECT_EQ(clipStart, s2.AlignedStart());    // queryStart (no soft clips left)
    EXPECT_EQ(clipEnd,   s2.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(tPos,      s2.ReferenceStart());  // 100 + startOffset
    EXPECT_EQ(tPos + 10, s2.ReferenceEnd());    // RefStart + 5M3D2M

    EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

    EXPECT_EQ(s2_seq_clipped,      s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_quals_clipped,    s2.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagBases_clipped, s2.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_tagQuals_clipped, s2.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_clipped, s2.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_clipped, s2.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_frames_clipped,   s2.IPD(Orientation::GENOMIC).Data());

    // s2 - REVERSE
    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
    EXPECT_EQ(clipStart, s2_rev.QueryStart());
    EXPECT_EQ(clipEnd,   s2_rev.QueryEnd());
    EXPECT_EQ(clipStart, s2_rev.AlignedStart());    // queryStart (no soft clips left)
    EXPECT_EQ(clipEnd,   s2_rev.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(tPos,      s2_rev.ReferenceStart());  // 100 + startOffset
    EXPECT_EQ(tPos + 10, s2_rev.ReferenceEnd());    // RefStart + 5M3D2M

    EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

    EXPECT_EQ(s2_seq_rev_clipped,      s2_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_quals_rev_clipped,    s2_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagBases_rev_clipped, s2_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_frames_rev_clipped,   s2_rev.IPD(Orientation::GENOMIC).Data());

    // s3 - FORWARD
    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(clipStart, s3.QueryStart());
    EXPECT_EQ(clipEnd,   s3.QueryEnd());
    EXPECT_EQ(clipStart, s3.AlignedStart());    // queryStart (no soft clips left)
    EXPECT_EQ(clipEnd,   s3.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(tPos,      s3.ReferenceStart());  // 100 + startOffset
    EXPECT_EQ(tPos + 8,  s3.ReferenceEnd());    // RefStart + 4M1D2D1M

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    EXPECT_EQ(s3_seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());

    // s3 - REVERSE
    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
    EXPECT_EQ(clipStart, s3_rev.QueryStart());
    EXPECT_EQ(clipEnd,   s3_rev.QueryEnd());
    EXPECT_EQ(clipStart, s3_rev.AlignedStart());    // queryStart (no soft clips left)
    EXPECT_EQ(clipEnd,   s3_rev.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(tPos,      s3_rev.ReferenceStart());  // 100 + startOffset
    EXPECT_EQ(tPos + 8,  s3_rev.ReferenceEnd());    // RefStart + 4M1D2D1M

    EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

    EXPECT_EQ(s3_seq_rev_clipped,      s3_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_rev_clipped,    s3_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_rev_clipped, s3_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_rev_clipped,   s3_rev.IPD(Orientation::GENOMIC).Data());
}

TEST(BamRecordClippingTest, ClipToReference_Basic)
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

    const string s1_cigar = "10M";
    const string s1_cigar_clipped = "5M";
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

    const string s2_cigar = "5M3D5M";
    const string s2_cigar_clipped = "3M2D";
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

    const string s3_cigar = "4M1D2I2D4M";
    const string s3_cigar_clipped = "2M1D2I2D";
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

    const BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
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

    // s0 - no clipping should have been done to unmapped record
    EXPECT_FALSE(s0.IsMapped());
    EXPECT_EQ(prototype.QueryStart(),     s0.QueryStart());
    EXPECT_EQ(prototype.QueryEnd(),       s0.QueryEnd());
    EXPECT_EQ(prototype.AlignedStart(),   s0.AlignedStart());
    EXPECT_EQ(prototype.AlignedEnd(),     s0.AlignedEnd());
    EXPECT_EQ(prototype.ReferenceStart(), s0.ReferenceStart());
    EXPECT_EQ(prototype.ReferenceEnd(),   s0.ReferenceEnd());
    EXPECT_EQ(prototype.Sequence(),       s0.Sequence());
    EXPECT_EQ(prototype.Qualities(),      s0.Qualities());
    EXPECT_EQ(prototype.DeletionTag(),    s0.DeletionTag());
    EXPECT_EQ(prototype.DeletionQV(),     s0.DeletionQV());
    EXPECT_EQ(prototype.LabelQV(),        s0.LabelQV());
    EXPECT_EQ(prototype.AltLabelQV(),     s0.AltLabelQV());
    EXPECT_EQ(prototype.IPD(),            s0.IPD());

    // s1 - FORWARD
    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
    EXPECT_EQ(502,   s1.QueryStart());
    EXPECT_EQ(507,   s1.QueryEnd());
    EXPECT_EQ(502,   s1.AlignedStart());       // queryStart (no soft clips)
    EXPECT_EQ(507,   s1.AlignedEnd());         // alignStart + seqLength
    EXPECT_EQ(clipStart, s1.ReferenceStart()); // clipStart
    EXPECT_EQ(clipEnd,   s1.ReferenceEnd());   // clipEnd

    EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

    EXPECT_EQ(s1_seq_clipped,      s1.Sequence());
    EXPECT_EQ(s1_quals_clipped,    s1.Qualities().Fastq());
    EXPECT_EQ(s1_tagBases_clipped, s1.DeletionTag());
    EXPECT_EQ(s1_tagQuals_clipped, s1.DeletionQV().Fastq());
    EXPECT_EQ(s1_tagQuals_clipped, s1.LabelQV().Fastq());
    EXPECT_EQ(s1_tagQuals_clipped, s1.AltLabelQV().Fastq());
    EXPECT_EQ(s1_frames_clipped,   s1.IPD().Data());

    // s1 - REVERSE
    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
    EXPECT_EQ(503, s1_rev.QueryStart());
    EXPECT_EQ(508, s1_rev.QueryEnd());
    EXPECT_EQ(503, s1_rev.AlignedStart());          // queryStart (no soft clips)
    EXPECT_EQ(508, s1_rev.AlignedEnd());            // alignStart + seqLength
    EXPECT_EQ(clipStart, s1_rev.ReferenceStart());  // clipStart
    EXPECT_EQ(clipEnd,   s1_rev.ReferenceEnd());    // clipEnd

    EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

    EXPECT_EQ(s1_seq_rev_clipped,      s1_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s1_quals_rev_clipped,    s1_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagBases_rev_clipped, s1_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_frames_rev_clipped,   s1_rev.IPD(Orientation::GENOMIC).Data());

    // s2 - FORWARD
    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
    EXPECT_EQ(502, s2.QueryStart());
    EXPECT_EQ(505, s2.QueryEnd());
    EXPECT_EQ(502, s2.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(505, s2.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(clipStart, s2.ReferenceStart());   // clipStart
    EXPECT_EQ(clipEnd,   s2.ReferenceEnd());     // clipEnd

    EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

    EXPECT_EQ(s2_seq_clipped,      s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_quals_clipped,    s2.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagBases_clipped, s2.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_tagQuals_clipped, s2.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_clipped, s2.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_clipped, s2.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_frames_clipped,   s2.IPD(Orientation::GENOMIC).Data());

    // s2 - REVERSE
    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
    EXPECT_EQ(505, s2_rev.QueryStart());
    EXPECT_EQ(508, s2_rev.QueryEnd());
    EXPECT_EQ(505, s2_rev.AlignedStart());    // queryStart (no soft clips)
    EXPECT_EQ(508, s2_rev.AlignedEnd());      // alignStart + seqLength
    EXPECT_EQ(clipStart, s2_rev.ReferenceStart());  // clipStart
    EXPECT_EQ(clipEnd,   s2_rev.ReferenceEnd());    // clipEnd

    EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

    EXPECT_EQ(s2_seq_rev_clipped,      s2_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_quals_rev_clipped,    s2_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagBases_rev_clipped, s2_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_frames_rev_clipped,   s2_rev.IPD(Orientation::GENOMIC).Data());

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

    EXPECT_EQ(s3_seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());

    // s3 - REVERSE
    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
    EXPECT_EQ(504, s3_rev.QueryStart());
    EXPECT_EQ(508, s3_rev.QueryEnd());
    EXPECT_EQ(504, s3_rev.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(508, s3_rev.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(clipStart, s3_rev.ReferenceStart());  // clipStart
    EXPECT_EQ(clipEnd,   s3_rev.ReferenceEnd());    // clipEnd

    EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

    EXPECT_EQ(s3_seq_rev_clipped,      s3_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_rev_clipped,    s3_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_rev_clipped, s3_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_rev_clipped,   s3_rev.IPD(Orientation::GENOMIC).Data());
}

TEST(BamRecordClippingTest, ClipToReference_WithSoftClips)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const string seq      = "TTAACCGTTAGCAAA";
    const string quals    = "--?]?]?]?]?*+++";
    const string tagBases = "TTAACCGTTAGCAAA";
    const string tagQuals = "--?]?]?]?]?*+++";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };

    const int32_t  tId     = 0;
    const Position tPos    = 100;
    const uint8_t  mapQual = 80;

    const Position clipStart = 102;
    const Position clipEnd   = 107;

    const string seq_rev      = "TTTGCTAACGGTTAA";
    const string quals_rev    = "+++*?]?]?]?]?--";
    const f_data frames_rev   = { 10, 10, 10, 20, 30, 10, 40, 40, 30, 20, 20, 10, 10, 40, 40 };

    const string s1_cigar = "2S10M3S";
    const string s1_cigar_clipped = "5M";
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

    const string s2_cigar = "2S5M3D5M3S";
    const string s2_cigar_clipped = "3M2D";
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

    const string s3_cigar = "2S4M1D2I2D4M3S";
    const string s3_cigar_clipped = "2M1D2I2D";
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

    const BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
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
    EXPECT_EQ(512,       s1.AlignedEnd());      // alignedStart + 10M
    EXPECT_EQ(tPos,      s1.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 10, s1.ReferenceEnd());    // tPos + 10M

    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(500,       s1_rev.QueryStart());      // queryStart
    EXPECT_EQ(515,       s1_rev.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(503,       s1_rev.AlignedStart());    // queryStart + 3S
    EXPECT_EQ(513,       s1_rev.AlignedEnd());      // alignedStart + 10M
    EXPECT_EQ(tPos,      s1_rev.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 10, s1_rev.ReferenceEnd());    // tPos + 10M

    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(500,       s2.QueryStart());      // queryStart
    EXPECT_EQ(515,       s2.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(502,       s2.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(512,       s2.AlignedEnd());      // alignedStart + 5M5M
    EXPECT_EQ(tPos,      s2.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 13, s2.ReferenceEnd());    // tPos + 5M3D5M

    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(500,       s2_rev.QueryStart());      // queryStart
    EXPECT_EQ(515,       s2_rev.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(503,       s2_rev.AlignedStart());    // queryStart + S
    EXPECT_EQ(513,       s2_rev.AlignedEnd());      // alignedStart + 5M5M
    EXPECT_EQ(tPos,      s2_rev.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 13, s2_rev.ReferenceEnd());    // tPos + 5M3D5M

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(500,       s3.QueryStart());      // queryStart
    EXPECT_EQ(515,       s3.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(502,       s3.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(512,       s3.AlignedEnd());      // alignedStart + 4M2I4M
    EXPECT_EQ(tPos,      s3.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 11, s3.ReferenceEnd());    // tPos + 4M1D2D4M

    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(500,       s3_rev.QueryStart());      // queryStart
    EXPECT_EQ(515,       s3_rev.QueryEnd());        // queryStart + seqLength
    EXPECT_EQ(503,       s3_rev.AlignedStart());    // queryStart + 2S
    EXPECT_EQ(513,       s3_rev.AlignedEnd());      // alignedStart + 4M2I4M
    EXPECT_EQ(tPos,      s3_rev.ReferenceStart());  // tPos
    EXPECT_EQ(tPos + 11, s3_rev.ReferenceEnd());    // tPos + 4M1D2D4M

    s0.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s1.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s2.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s3.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s1_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s2_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);
    s3_rev.Clip(ClipType::CLIP_TO_REFERENCE, clipStart, clipEnd);

    // s0 - no clipping should have been done to unmapped record
    EXPECT_FALSE(s0.IsMapped());
    EXPECT_EQ(prototype.QueryStart(),     s0.QueryStart());
    EXPECT_EQ(prototype.QueryEnd(),       s0.QueryEnd());
    EXPECT_EQ(prototype.AlignedStart(),   s0.AlignedStart());
    EXPECT_EQ(prototype.AlignedEnd(),     s0.AlignedEnd());
    EXPECT_EQ(prototype.ReferenceStart(), s0.ReferenceStart());
    EXPECT_EQ(prototype.ReferenceEnd(),   s0.ReferenceEnd());
    EXPECT_EQ(prototype.Sequence(),       s0.Sequence());
    EXPECT_EQ(prototype.Qualities(),      s0.Qualities());
    EXPECT_EQ(prototype.DeletionTag(),    s0.DeletionTag());
    EXPECT_EQ(prototype.DeletionQV(),     s0.DeletionQV());
    EXPECT_EQ(prototype.LabelQV(),        s0.LabelQV());
    EXPECT_EQ(prototype.AltLabelQV(),     s0.AltLabelQV());
    EXPECT_EQ(prototype.IPD(),            s0.IPD());

    // s1 - FORWARD
    EXPECT_TRUE(s1.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
    EXPECT_EQ(504,   s1.QueryStart());         // new queryStart
    EXPECT_EQ(509,   s1.QueryEnd());           // queryStart + new seqLength
    EXPECT_EQ(504,   s1.AlignedStart());       // queryStart (no soft clips remaining)
    EXPECT_EQ(509,   s1.AlignedEnd());         // alignStart + new seqLength
    EXPECT_EQ(clipStart, s1.ReferenceStart()); // clipStart
    EXPECT_EQ(clipEnd,   s1.ReferenceEnd());   // clipEnd

    EXPECT_EQ(s1_cigar_clipped, s1.CigarData().ToStdString());

    EXPECT_EQ(s1_seq_clipped,      s1.Sequence());
    EXPECT_EQ(s1_quals_clipped,    s1.Qualities().Fastq());
    EXPECT_EQ(s1_tagBases_clipped, s1.DeletionTag());
    EXPECT_EQ(s1_tagQuals_clipped, s1.DeletionQV().Fastq());
    EXPECT_EQ(s1_tagQuals_clipped, s1.LabelQV().Fastq());
    EXPECT_EQ(s1_tagQuals_clipped, s1.AltLabelQV().Fastq());
    EXPECT_EQ(s1_frames_clipped,   s1.IPD().Data());

    // s1 - REVERSE
    EXPECT_TRUE(s1_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
    EXPECT_EQ(506,   s1_rev.QueryStart());         // new queryStart
    EXPECT_EQ(511,   s1_rev.QueryEnd());           // queryStart + new seqLength
    EXPECT_EQ(506,   s1_rev.AlignedStart());       // queryStart (no soft clips remaining)
    EXPECT_EQ(511,   s1_rev.AlignedEnd());         // alignStart + new seqLength
    EXPECT_EQ(clipStart, s1_rev.ReferenceStart()); // clipStart
    EXPECT_EQ(clipEnd,   s1_rev.ReferenceEnd());   // clipEnd

    EXPECT_EQ(s1_cigar_clipped, s1_rev.CigarData().ToStdString());

    EXPECT_EQ(s1_seq_rev_clipped,      s1_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s1_quals_rev_clipped,    s1_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagBases_rev_clipped, s1_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_tagQuals_rev_clipped, s1_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s1_frames_rev_clipped,   s1_rev.IPD(Orientation::GENOMIC).Data());

    // s2 - FORWARD
    EXPECT_TRUE(s2.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
    EXPECT_EQ(504, s2.QueryStart());
    EXPECT_EQ(507, s2.QueryEnd());
    EXPECT_EQ(504, s2.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(507, s2.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(clipStart, s2.ReferenceStart());   // clipStart
    EXPECT_EQ(clipEnd,   s2.ReferenceEnd());     // clipEnd

    EXPECT_EQ(s2_cigar_clipped, s2.CigarData().ToStdString());

    EXPECT_EQ(s2_seq_clipped,      s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_quals_clipped,    s2.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagBases_clipped, s2.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_tagQuals_clipped, s2.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_clipped, s2.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_clipped, s2.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_frames_clipped,   s2.IPD(Orientation::GENOMIC).Data());

    // s2 - REVERSE
    EXPECT_TRUE(s2_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
    EXPECT_EQ(508,   s2_rev.QueryStart());         // new queryStart
    EXPECT_EQ(511,   s2_rev.QueryEnd());           // queryStart + new seqLength
    EXPECT_EQ(508,   s2_rev.AlignedStart());       // queryStart (no soft clips remaining)
    EXPECT_EQ(511,   s2_rev.AlignedEnd());         // alignStart + new seqLength
    EXPECT_EQ(clipStart, s2_rev.ReferenceStart()); // clipStart
    EXPECT_EQ(clipEnd,   s2_rev.ReferenceEnd());   // clipEnd

    EXPECT_EQ(s2_cigar_clipped, s2_rev.CigarData().ToStdString());

    EXPECT_EQ(s2_seq_rev_clipped,      s2_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_quals_rev_clipped,    s2_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagBases_rev_clipped, s2_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_tagQuals_rev_clipped, s2_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_frames_rev_clipped,   s2_rev.IPD(Orientation::GENOMIC).Data());

    // s3 - FORWARD
    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(504, s3.QueryStart());
    EXPECT_EQ(508, s3.QueryEnd());
    EXPECT_EQ(504, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(508, s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(clipStart, s3.ReferenceStart());   // clipStart
    EXPECT_EQ(clipEnd,   s3.ReferenceEnd());     // clipEnd

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    EXPECT_EQ(s3_seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());

    // s3 - REVERSE
    EXPECT_TRUE(s3_rev.IsMapped());
    EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
    EXPECT_EQ(507,   s3_rev.QueryStart());         // new queryStart
    EXPECT_EQ(511,   s3_rev.QueryEnd());           // queryStart + new seqLength
    EXPECT_EQ(507,   s3_rev.AlignedStart());       // queryStart (no soft clips remaining)
    EXPECT_EQ(511,   s3_rev.AlignedEnd());         // alignStart + new seqLength
    EXPECT_EQ(clipStart, s3_rev.ReferenceStart()); // clipStart
    EXPECT_EQ(clipEnd,   s3_rev.ReferenceEnd());   // clipEnd

    EXPECT_EQ(s3_cigar_clipped, s3_rev.CigarData().ToStdString());

    EXPECT_EQ(s3_seq_rev_clipped,      s3_rev.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_rev_clipped,    s3_rev.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_rev_clipped, s3_rev.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_rev_clipped, s3_rev.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_rev_clipped,   s3_rev.IPD(Orientation::GENOMIC).Data());
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

    const string s3_cigar = "4M1D2I2D4M";
    const string s3_cigar_clipped = "2M1D2I2D3M";

    BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    BamRecord s3 = prototype.Clipped(ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(clipStart, s3.QueryStart());
    EXPECT_EQ(clipEnd,   s3.QueryEnd());
    EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5M + 3D

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    EXPECT_EQ(seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());
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

    const string s3_cigar = "4M1D2I2D4M";
    const string s3_cigar_clipped = "2M1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };

    BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
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

    EXPECT_EQ(s3_seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());
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

    const string s3_cigar = "4M1D2I2D4M";
    const string s3_cigar_clipped = "2M1D2I2D3M";

    BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    prototype.Map(tId, tPos, Strand::FORWARD, s3_cigar, mapQual);

    const BamRecord s3 = BamRecord::Clipped(prototype, ClipType::CLIP_TO_QUERY, clipStart, clipEnd);

    EXPECT_TRUE(s3.IsMapped());
    EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
    EXPECT_EQ(clipStart, s3.QueryStart());
    EXPECT_EQ(clipEnd,   s3.QueryEnd());
    EXPECT_EQ(clipStart, s3.AlignedStart());     // queryStart (no soft clips)
    EXPECT_EQ(clipEnd,   s3.AlignedEnd());       // alignStart + seqLength
    EXPECT_EQ(102, s3.ReferenceStart());         // 100 + startOffset
    EXPECT_EQ(110, s3.ReferenceEnd());           // RefStart + 5M + 3D

    EXPECT_EQ(s3_cigar_clipped, s3.CigarData().ToStdString());

    EXPECT_EQ(seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());
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

    const string s3_cigar = "4M1D2I2D4M";
    const string s3_cigar_clipped = "2M1D2I2D";
    const string s3_seq_clipped      = "CCGT";
    const string s3_quals_clipped    = "?]?]";
    const string s3_tagBases_clipped = s3_seq_clipped;
    const string s3_tagQuals_clipped = s3_quals_clipped;
    const f_data s3_frames_clipped   = { 20, 20, 30, 40 };

    BamRecord prototype = tests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
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

    EXPECT_EQ(s3_seq_clipped,      s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_quals_clipped,    s3.Qualities(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagBases_clipped, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_tagQuals_clipped, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.LabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_tagQuals_clipped, s3.AltLabelQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_frames_clipped,   s3.IPD(Orientation::GENOMIC).Data());
}

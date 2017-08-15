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

#include <chrono>
#include <string>

#include <gtest/gtest.h>

#ifdef PBBAM_TESTING
#define private public
#endif

#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordView.h>
#include <pbbam/BamTagCodec.h>

using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

typedef vector<uint16_t> f_data;

namespace BamRecordMappingTests {

static
BamRecord MakeRecord(const Position qStart,
                     const Position qEnd,
                     const string& seq,
                     const string& quals,
                     const string& tagBases,
                     const string& tagQuals,
                     const f_data& frames)
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

} // namespace BamRecordMappingTests

TEST(BamRecordMappingTest, BasicMap)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const uint8_t mapQual = 80;

    const string seq_rev   = "GCTAACGGTT";
    const string quals_rev = "*?]?]?]?]?";
    const string tagBases_rev = seq_rev;
    const string tagQuals_rev = quals_rev;
    const f_data frames_rev = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const string s1_cigar = "10=";
    const string s2_cigar = "5=3D5=";
    const string s3_cigar = "4=1D2I2D4=";

    BamRecord s1 = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s2 = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s3 = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s1_rev = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s2_rev = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s3_rev = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);

    s1.Map(0, 100, Strand::FORWARD, s1_cigar, mapQual);
    s2.Map(0, 100, Strand::FORWARD, s2_cigar, mapQual);
    s3.Map(0, 100, Strand::FORWARD, s3_cigar, mapQual);
    s1_rev.Map(0, 100, Strand::REVERSE, s1_cigar, mapQual);
    s2_rev.Map(0, 100, Strand::REVERSE, s2_cigar, mapQual);
    s3_rev.Map(0, 100, Strand::REVERSE, s3_cigar, mapQual);

    {   // s1 - FORWARD
        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(0, s1.ReferenceId());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(mapQual, s1.MapQuality());

        EXPECT_EQ(qStart, s1.QueryStart());
        EXPECT_EQ(qEnd,   s1.QueryEnd());
        EXPECT_EQ(500, s1.AlignedStart());
        EXPECT_EQ(510, s1.AlignedEnd());         // 500 + 10=
        EXPECT_EQ(100, s1.ReferenceStart());
        EXPECT_EQ(110, s1.ReferenceEnd());       // 100 + 10=

        const BamRecordView view
        {
            s1,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq,      view.Sequence());
        EXPECT_EQ(quals,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases, view.DeletionTags());
        EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   view.IPD().Data());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(0, s1_rev.ReferenceId());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(mapQual, s1_rev.MapQuality());

        EXPECT_EQ(qStart, s1_rev.QueryStart());
        EXPECT_EQ(qEnd,   s1_rev.QueryEnd());
        EXPECT_EQ(500, s1_rev.AlignedStart());
        EXPECT_EQ(510, s1_rev.AlignedEnd());         // 500 + 10=
        EXPECT_EQ(100, s1_rev.ReferenceStart());
        EXPECT_EQ(110, s1_rev.ReferenceEnd());       // 100 + 10=

        // native
        const BamRecordView nativeView
        {
            s1_rev,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq,      nativeView.Sequence());
        EXPECT_EQ(quals,    nativeView.Qualities().Fastq());
        EXPECT_EQ(tagBases, nativeView.DeletionTags());
        EXPECT_EQ(tagQuals, nativeView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   nativeView.IPD().Data());

        // - genomic
        const BamRecordView genomicView
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq_rev,      genomicView.Sequence());
        EXPECT_EQ(quals_rev,    genomicView.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev, genomicView.DeletionTags());
        EXPECT_EQ(tagQuals_rev, genomicView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev,   genomicView.IPD().Data());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(0, s2.ReferenceId());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(mapQual, s2.MapQuality());

        EXPECT_EQ(qStart, s2.QueryStart());
        EXPECT_EQ(qEnd,   s2.QueryEnd());
        EXPECT_EQ(500, s2.AlignedStart());
        EXPECT_EQ(510, s2.AlignedEnd());         // 500 + 10=
        EXPECT_EQ(100, s2.ReferenceStart());
        EXPECT_EQ(113, s2.ReferenceEnd());      // 100 + 10= + 3D

        const BamRecordView view
        {
            s2,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq,      view.Sequence());
        EXPECT_EQ(quals,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases, view.DeletionTags());
        EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(0, s2_rev.ReferenceId());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(mapQual, s2_rev.MapQuality());

        EXPECT_EQ(qStart, s2_rev.QueryStart());
        EXPECT_EQ(qEnd,   s2_rev.QueryEnd());
        EXPECT_EQ(500, s2_rev.AlignedStart());
        EXPECT_EQ(510, s2_rev.AlignedEnd());         // 500 + 10=
        EXPECT_EQ(100, s2_rev.ReferenceStart());
        EXPECT_EQ(113, s2_rev.ReferenceEnd());      // 100 + 10= + 3D

        // - native
        const BamRecordView nativeView
        {
            s2_rev,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq,      nativeView.Sequence());
        EXPECT_EQ(quals,    nativeView.Qualities().Fastq());
        EXPECT_EQ(tagBases, nativeView.DeletionTags());
        EXPECT_EQ(tagQuals, nativeView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   nativeView.IPD().Data());

        // - genomic
        const BamRecordView genomicView
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq_rev,      genomicView.Sequence());
        EXPECT_EQ(quals_rev,    genomicView.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev, genomicView.DeletionTags());
        EXPECT_EQ(tagQuals_rev, genomicView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev,   genomicView.IPD().Data());
    }

    {   // s3 - FORWARD

        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(0, s3.ReferenceId());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(mapQual, s3.MapQuality());

        EXPECT_EQ(qStart, s3.QueryStart());
        EXPECT_EQ(qEnd,   s3.QueryEnd());
        EXPECT_EQ(500, s3.AlignedStart());
        EXPECT_EQ(510, s3.AlignedEnd());         // 500 + 8= + 2I
        EXPECT_EQ(100, s3.ReferenceStart());
        EXPECT_EQ(111, s3.ReferenceEnd());      // 100 + 8= + 3D

        const BamRecordView view
        {
            s3,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq,      view.Sequence());
        EXPECT_EQ(quals,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases, view.DeletionTags());
        EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   view.IPD().Data());
    }

    {   // s3 - REVERSE

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(0, s3_rev.ReferenceId());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(mapQual, s3_rev.MapQuality());

        EXPECT_EQ(qStart, s3_rev.QueryStart());
        EXPECT_EQ(qEnd,   s3_rev.QueryEnd());
        EXPECT_EQ(500, s3_rev.AlignedStart());
        EXPECT_EQ(510, s3_rev.AlignedEnd());         // 500 + 8= + 2I
        EXPECT_EQ(100, s3_rev.ReferenceStart());
        EXPECT_EQ(111, s3_rev.ReferenceEnd());      // 100 + 8= + 3D

        // - native
        const BamRecordView nativeView
        {
            s3_rev,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq,      nativeView.Sequence());
        EXPECT_EQ(quals,    nativeView.Qualities().Fastq());
        EXPECT_EQ(tagBases, nativeView.DeletionTags());
        EXPECT_EQ(tagQuals, nativeView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   nativeView.IPD().Data());

        // - genomic
        const BamRecordView genomicView
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq_rev,      genomicView.Sequence());
        EXPECT_EQ(quals_rev,    genomicView.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev, genomicView.DeletionTags());
        EXPECT_EQ(tagQuals_rev, genomicView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev,   genomicView.IPD().Data());
    }
}

TEST(BamRecordMappingTest, SoftClipMapping)
{
    const Position qStart = 500;
    const Position qEnd   = 515;
    const string seq      = "TTAACCGTTAGCAAA";
    const string quals    = "--?]?]?]?]?*+++";
    const string tagBases = "TTAACCGTTAGCAAA";
    const string tagQuals = "--?]?]?]?]?*+++";
    const f_data frames   = { 40, 40, 10, 10, 20, 20, 30, 40, 40, 10, 30, 20, 10, 10, 10 };
    const uint8_t mapQual = 80;

    const string clipped_seq   = "AACCGTTAGC";
    const string clipped_quals = "?]?]?]?]?*";
    const string clipped_tagBases   = "AACCGTTAGC";
    const string clipped_tagQuals = "?]?]?]?]?*";
    const f_data clipped_frames = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };

    const string seq_rev   = "TTTGCTAACGGTTAA";
    const string quals_rev = "+++*?]?]?]?]?--";
    const string tagBases_rev = seq_rev;
    const string tagQuals_rev = quals_rev;
    const f_data frames_rev = { 10, 10, 10, 20, 30, 10, 40, 40, 30, 20, 20, 10, 10, 40, 40 };

    const string clipped_seq_rev   = "GCTAACGGTT";
    const string clipped_quals_rev = "*?]?]?]?]?";
    const string clipped_tagBases_rev = clipped_seq_rev;
    const string clipped_tagQuals_rev = clipped_quals_rev;
    const f_data clipped_frames_rev = { 20, 30, 10, 40, 40, 30, 20, 20, 10, 10 };

    const string s1_cigar = "2S10=3S";
    const string s2_cigar = "2S5=3D5=3S";
    const string s3_cigar = "2S4=1D2I2D4=3S";

    BamRecord s1 = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s2 = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s3 = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s1_rev = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s2_rev = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    BamRecord s3_rev = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);

    s1.Map(0, 100, Strand::FORWARD, s1_cigar, mapQual);
    s2.Map(0, 100, Strand::FORWARD, s2_cigar, mapQual);
    s3.Map(0, 100, Strand::FORWARD, s3_cigar, mapQual);
    s1_rev.Map(0, 100, Strand::REVERSE, s1_cigar, mapQual);
    s2_rev.Map(0, 100, Strand::REVERSE, s2_cigar, mapQual);
    s3_rev.Map(0, 100, Strand::REVERSE, s3_cigar, mapQual);

    {   // s1 - FORWARD

        EXPECT_TRUE(s1.IsMapped());
        EXPECT_EQ(0, s1.ReferenceId());
        EXPECT_EQ(Strand::FORWARD, s1.AlignedStrand());
        EXPECT_EQ(mapQual, s1.MapQuality());

        EXPECT_EQ(qStart, s1.QueryStart());      // 500
        EXPECT_EQ(qEnd,   s1.QueryEnd());        // QStart + seqLength
        EXPECT_EQ(502, s1.AlignedStart());       // QStart + 2S
        EXPECT_EQ(512, s1.AlignedEnd());         // AStart + 10=
        EXPECT_EQ(100, s1.ReferenceStart());     // 100
        EXPECT_EQ(110, s1.ReferenceEnd());       // RefStart + 10=

        const BamRecordView view
        {
            s1,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq,      view.Sequence());
        EXPECT_EQ(quals,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases, view.DeletionTags());
        EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   view.IPD().Data());
    }

    {   // s1 - REVERSE

        EXPECT_TRUE(s1_rev.IsMapped());
        EXPECT_EQ(0, s1_rev.ReferenceId());
        EXPECT_EQ(Strand::REVERSE, s1_rev.AlignedStrand());
        EXPECT_EQ(mapQual, s1_rev.MapQuality());

        EXPECT_EQ(qStart, s1_rev.QueryStart());      // 500
        EXPECT_EQ(qEnd,   s1_rev.QueryEnd());        // QStart + seqLength
        EXPECT_EQ(503, s1_rev.AlignedStart());       // QStart + 3S
        EXPECT_EQ(513, s1_rev.AlignedEnd());         // AStart + 10=
        EXPECT_EQ(100, s1_rev.ReferenceStart());     // 100
        EXPECT_EQ(110, s1_rev.ReferenceEnd());       // RefStart + 10=

        // - native
        const BamRecordView nativeView
        {
            s1_rev,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq,      nativeView.Sequence());
        EXPECT_EQ(quals,    nativeView.Qualities().Fastq());
        EXPECT_EQ(tagBases, nativeView.DeletionTags());
        EXPECT_EQ(tagQuals, nativeView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   nativeView.IPD().Data());

        // - genomic
        const BamRecordView genomicView
        {
            s1_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq_rev,      genomicView.Sequence());
        EXPECT_EQ(quals_rev,    genomicView.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev, genomicView.DeletionTags());
        EXPECT_EQ(tagQuals_rev, genomicView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev,   genomicView.IPD().Data());
    }

    {   // s2 - FORWARD

        EXPECT_TRUE(s2.IsMapped());
        EXPECT_EQ(0, s2.ReferenceId());
        EXPECT_EQ(Strand::FORWARD, s2.AlignedStrand());
        EXPECT_EQ(mapQual, s2.MapQuality());

        EXPECT_EQ(qStart, s2.QueryStart());      // 500
        EXPECT_EQ(qEnd,   s2.QueryEnd());        // QStart + seqLength
        EXPECT_EQ(502, s2.AlignedStart());       // QStart + 2S
        EXPECT_EQ(512, s2.AlignedEnd());         // AStart + 10=
        EXPECT_EQ(100, s2.ReferenceStart());     // 100
        EXPECT_EQ(113, s2.ReferenceEnd());       // RefStart + 10= + 3D

        const BamRecordView view
        {
            s2,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq,      view.Sequence());
        EXPECT_EQ(quals,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases, view.DeletionTags());
        EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   view.IPD().Data());
    }

    {   // s2 - REVERSE

        EXPECT_TRUE(s2_rev.IsMapped());
        EXPECT_EQ(0, s2_rev.ReferenceId());
        EXPECT_EQ(Strand::REVERSE, s2_rev.AlignedStrand());
        EXPECT_EQ(mapQual, s2_rev.MapQuality());

        EXPECT_EQ(qStart, s2_rev.QueryStart());      // 500
        EXPECT_EQ(qEnd,   s2_rev.QueryEnd());        // QStart + seqLength
        EXPECT_EQ(503, s2_rev.AlignedStart());       // QStart + 3S
        EXPECT_EQ(513, s2_rev.AlignedEnd());         // AStart + 10=
        EXPECT_EQ(100, s2_rev.ReferenceStart());     // 100
        EXPECT_EQ(113, s2_rev.ReferenceEnd());       // RefStart + 10= + 3D

        // - native
        const BamRecordView nativeView
        {
            s2_rev,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq,      nativeView.Sequence());
        EXPECT_EQ(quals,    nativeView.Qualities().Fastq());
        EXPECT_EQ(tagBases, nativeView.DeletionTags());
        EXPECT_EQ(tagQuals, nativeView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   nativeView.IPD().Data());

        // - genomic
        const BamRecordView genomicView
        {
            s2_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq_rev,      genomicView.Sequence());
        EXPECT_EQ(quals_rev,    genomicView.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev, genomicView.DeletionTags());
        EXPECT_EQ(tagQuals_rev, genomicView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev,   genomicView.IPD().Data());
    }

    {   // s3 - FORWARD

        EXPECT_TRUE(s3.IsMapped());
        EXPECT_EQ(0, s3.ReferenceId());
        EXPECT_EQ(Strand::FORWARD, s3.AlignedStrand());
        EXPECT_EQ(mapQual, s3.MapQuality());

        EXPECT_EQ(qStart, s3.QueryStart());      // 500
        EXPECT_EQ(qEnd,   s3.QueryEnd());        // QStart + seqLength
        EXPECT_EQ(502, s3.AlignedStart());       // QStart + 2S
        EXPECT_EQ(512, s3.AlignedEnd());         // AStart + 8= + 2I
        EXPECT_EQ(100, s3.ReferenceStart());     // 100
        EXPECT_EQ(111, s3.ReferenceEnd());       // RefStart + 8= + 3D

        const BamRecordView view
        {
            s2,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };

        EXPECT_EQ(seq,      view.Sequence());
        EXPECT_EQ(quals,    view.Qualities().Fastq());
        EXPECT_EQ(tagBases, view.DeletionTags());
        EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   view.IPD().Data());
    }

    {   // s3 - REVERSE

        EXPECT_TRUE(s3_rev.IsMapped());
        EXPECT_EQ(0, s3_rev.ReferenceId());
        EXPECT_EQ(Strand::REVERSE, s3_rev.AlignedStrand());
        EXPECT_EQ(mapQual, s3_rev.MapQuality());

        EXPECT_EQ(qStart, s3_rev.QueryStart());      // 500
        EXPECT_EQ(qEnd,   s3_rev.QueryEnd());        // QStart + seqLength
        EXPECT_EQ(503, s3_rev.AlignedStart());       // QStart + 3S
        EXPECT_EQ(513, s3_rev.AlignedEnd());         // AStart + 8= + 2I
        EXPECT_EQ(100, s3_rev.ReferenceStart());     // 100
        EXPECT_EQ(111, s3_rev.ReferenceEnd());       // RefStart + 8= + 3D

        // - native
        const BamRecordView nativeView
        {
            s3_rev,
            Orientation::NATIVE,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq,      nativeView.Sequence());
        EXPECT_EQ(quals,    nativeView.Qualities().Fastq());
        EXPECT_EQ(tagBases, nativeView.DeletionTags());
        EXPECT_EQ(tagQuals, nativeView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals, nativeView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames,   nativeView.IPD().Data());

        // - genomic
        const BamRecordView genomicView
        {
            s3_rev,
            Orientation::GENOMIC,
            false,
            false,
            PulseBehavior::ALL
        };
        EXPECT_EQ(seq_rev,      genomicView.Sequence());
        EXPECT_EQ(quals_rev,    genomicView.Qualities().Fastq());
        EXPECT_EQ(tagBases_rev, genomicView.DeletionTags());
        EXPECT_EQ(tagQuals_rev, genomicView.DeletionQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.LabelQVs().Fastq());
        EXPECT_EQ(tagQuals_rev, genomicView.AltLabelQVs().Fastq());
        EXPECT_EQ(frames_rev,   genomicView.IPD().Data());
    }
}

TEST(BamRecordMappingTest, MappedCopy)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const uint8_t mapQual = 80;
    const string cigar    = "4=1D2I2D4=";

    const BamRecord orig = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    const BamRecord mapped = orig.Mapped(0, 100, Strand::FORWARD, cigar, mapQual);

    EXPECT_TRUE(mapped.IsMapped());
    EXPECT_EQ(0, mapped.ReferenceId());
    EXPECT_EQ(Strand::FORWARD, mapped.AlignedStrand());
    EXPECT_EQ(mapQual, mapped.MapQuality());

    EXPECT_EQ(500, mapped.QueryStart());      // 500
    EXPECT_EQ(510, mapped.QueryEnd());        // QStart + seqLength
    EXPECT_EQ(500, mapped.AlignedStart());    // QStart
    EXPECT_EQ(510, mapped.AlignedEnd());      // QStart + 8= + 2I
    EXPECT_EQ(100, mapped.ReferenceStart());  // 100
    EXPECT_EQ(111, mapped.ReferenceEnd());    // RefStart + 8= + 3D

    const BamRecordView view
    {
        mapped,
        Orientation::NATIVE,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(seq,      view.Sequence());
    EXPECT_EQ(quals,    view.Qualities().Fastq());
    EXPECT_EQ(tagBases, view.DeletionTags());
    EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
    EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
    EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
    EXPECT_EQ(frames,   view.IPD().Data());
}

TEST(BamRecordMappingTest, StaticMapped)
{
    const Position qStart = 500;
    const Position qEnd   = 510;
    const string seq      = "AACCGTTAGC";
    const string quals    = "?]?]?]?]?*";
    const string tagBases = "AACCGTTAGC";
    const string tagQuals = "?]?]?]?]?*";
    const f_data frames   = { 10, 10, 20, 20, 30, 40, 40, 10, 30, 20 };
    const uint8_t mapQual = 80;
    const string cigar    = "4=1D2I2D4=";

    const BamRecord orig = BamRecordMappingTests::MakeRecord(qStart, qEnd, seq, quals, tagBases, tagQuals, frames);
    const BamRecord mapped = BamRecord::Mapped(orig, 0, 100, Strand::FORWARD, cigar, mapQual);

    EXPECT_TRUE(mapped.IsMapped());
    EXPECT_EQ(0, mapped.ReferenceId());
    EXPECT_EQ(Strand::FORWARD, mapped.AlignedStrand());
    EXPECT_EQ(mapQual, mapped.MapQuality());

    EXPECT_EQ(500, mapped.QueryStart());      // 500
    EXPECT_EQ(510, mapped.QueryEnd());        // QStart + seqLength
    EXPECT_EQ(500, mapped.AlignedStart());    // QStart
    EXPECT_EQ(510, mapped.AlignedEnd());      // QStart + 8= + 2I
    EXPECT_EQ(100, mapped.ReferenceStart());  // 100
    EXPECT_EQ(111, mapped.ReferenceEnd());    // RefStart + 8= + 3D

    const BamRecordView view
    {
        mapped,
        Orientation::NATIVE,
        false,
        false,
        PulseBehavior::ALL
    };

    EXPECT_EQ(seq,      view.Sequence());
    EXPECT_EQ(quals,    view.Qualities().Fastq());
    EXPECT_EQ(tagBases, view.DeletionTags());
    EXPECT_EQ(tagQuals, view.DeletionQVs().Fastq());
    EXPECT_EQ(tagQuals, view.LabelQVs().Fastq());
    EXPECT_EQ(tagQuals, view.AltLabelQVs().Fastq());
    EXPECT_EQ(frames,   view.IPD().Data());
}

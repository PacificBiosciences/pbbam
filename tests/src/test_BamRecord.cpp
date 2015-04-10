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
#include <string>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace tests {

static
BamRecordImpl CreateBamImpl(void)
{
    BamRecordImpl bam;
    bam.Bin(42);
    bam.Flag(42);
    bam.InsertSize(42);
    bam.MapQuality(42);
    bam.MatePosition(42);
    bam.MateReferenceId(42);
    bam.Position(42);
    bam.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam.Tags(tags);

    return bam;
}

static
BamRecord CreateBam(void)
{
    BamRecord bam;
    bam.impl_ = CreateBamImpl();
    return bam;
}

static
void CheckRawData(const BamRecordImpl& bam)
{
    // ensure raw data (lengths at least) matches API-facing data

    const uint32_t expectedNameLength  = bam.Name().size() + 1;
    const uint32_t expectedNumCigarOps = bam.CigarData().size();
    const int32_t  expectedSeqLength   = bam.Sequence().length();
    const size_t   expectedTagsLength  = BamTagCodec::Encode(bam.Tags()).size();

    //  Name        CIGAR         Sequence       Quals      Tags
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + << TAGS >>

    const int expectedTotalDataLength = expectedNameLength +
                                        (expectedNumCigarOps * 4) +
                                        (expectedSeqLength+1)/2 +
                                         expectedSeqLength +
                                         expectedTagsLength;

    EXPECT_TRUE((bool)bam.d_);
    EXPECT_EQ(expectedNameLength,      bam.d_->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps,     bam.d_->core.n_cigar);
    EXPECT_EQ(expectedSeqLength,       bam.d_->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, bam.d_->l_data);
}

static
void CheckRawData(const BamRecord& bam)
{ CheckRawData(bam.impl_); }

static
BamRecordImpl MakeCigaredImpl(const string& seq,
                              const string& cigar,
                              const bool isReverseStrand)
{
    BamRecordImpl impl;
    impl.SetMapped(true).ReferenceId(0).Position(0).MapQuality(0);
    impl.CigarData( Cigar::FromStdString(cigar) );
    impl.MateReferenceId(-1).MatePosition(-1).InsertSize(0);
    impl.SetSequenceAndQualities(seq, string(seq.size(), '*'));
    impl.SetReverseStrand(isReverseStrand);
    return impl;
}

static
BamRecord MakeCigaredRecord(const string& seq,
                            const string& cigar,
                            const bool isReverseStrand)
{
    const BamRecordImpl impl = MakeCigaredImpl(seq, cigar, isReverseStrand);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredBaseRecord(const string& bases,
                                const string& cigar,
                                const bool isReverseStrand)
{
    const string seq = string(bases.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, isReverseStrand);

    TagCollection tags;
    tags["dt"] = bases;
    tags["st"] = bases;
    impl.Tags(tags);

    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredFrameRecord(const vector<uint16_t>& frames,
                                 const string& cigar,
                                 const bool isReverseStrand)
{
    const string seq = string(frames.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, isReverseStrand);

    TagCollection tags;
    tags["ip"] = frames;
    tags["pw"] = frames;
    impl.Tags(tags);

    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredQualRecord(const string& quals,
                                const string& cigar,
                                const bool isReverseStrand)
{
    const string seq = string(quals.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, isReverseStrand);

    TagCollection tags;
    tags["dq"] = quals;
    tags["iq"] = quals;
    tags["mq"] = quals;
    tags["sq"] = quals;
    impl.Tags(tags);

    return BamRecord(std::move(impl));
}

} // namespace tests

TEST(BamRecordTest, DefaultValues)
{
    BamRecord bam;
    const string emptyString;

    // BamRecordImpl data
    EXPECT_EQ(0, bam.impl_.Bin());
    EXPECT_EQ(BamRecordImpl::UNMAPPED, bam.impl_.Flag());  // forced init unmapped
    EXPECT_EQ(0, bam.impl_.InsertSize());
    EXPECT_EQ(0, bam.impl_.MapQuality());
    EXPECT_EQ(-1, bam.impl_.MateReferenceId());
    EXPECT_EQ(-1, bam.impl_.MatePosition());
    EXPECT_EQ(-1, bam.impl_.Position());
    EXPECT_EQ(-1, bam.impl_.ReferenceId());
    EXPECT_EQ(0, bam.impl_.Tags().size());

    EXPECT_FALSE(bam.impl_.IsDuplicate());
    EXPECT_FALSE(bam.impl_.IsFailedQC());
    EXPECT_FALSE(bam.impl_.IsFirstMate());
    EXPECT_FALSE(bam.impl_.IsMapped());             // forced init unmapped
    EXPECT_TRUE(bam.impl_.IsMateMapped());
    EXPECT_FALSE(bam.impl_.IsMateReverseStrand());
    EXPECT_FALSE(bam.impl_.IsPaired());
    EXPECT_TRUE(bam.impl_.IsPrimaryAlignment());
    EXPECT_FALSE(bam.impl_.IsProperPair());
    EXPECT_FALSE(bam.impl_.IsReverseStrand());
    EXPECT_FALSE(bam.impl_.IsSecondMate());
    EXPECT_FALSE(bam.impl_.IsSupplementaryAlignment());

    EXPECT_EQ(emptyString, bam.impl_.Name());
    EXPECT_EQ(emptyString, bam.impl_.CigarData().ToStdString());
    EXPECT_EQ(emptyString, bam.impl_.Sequence());
    EXPECT_EQ(emptyString, bam.impl_.Qualities().Fastq());

    // PacBio data
    EXPECT_EQ(-1, bam.AlignedStart());
    EXPECT_EQ(-1, bam.AlignedEnd());
    EXPECT_EQ(-1, bam.HoleNumber());
    EXPECT_EQ(-1, bam.NumPasses());
    EXPECT_EQ(-1, bam.QueryStart());
    EXPECT_EQ(-1, bam.QueryEnd());
    EXPECT_EQ(0, bam.ReadAccuracy());
    EXPECT_EQ(-1, bam.ReferenceStart());
    EXPECT_EQ(-1, bam.ReferenceEnd());

    EXPECT_TRUE(bam.DeletionQV().empty());
    EXPECT_TRUE(bam.DeletionTag().empty());
    EXPECT_TRUE(bam.InsertionQV().empty());
    EXPECT_TRUE(bam.MergeQV().empty());
    EXPECT_TRUE(bam.SubstitutionQV().empty());
    EXPECT_TRUE(bam.SubstitutionTag().empty());

    EXPECT_FALSE(bam.HasDeletionQV());
    EXPECT_FALSE(bam.HasDeletionTag());
    EXPECT_FALSE(bam.HasInsertionQV());
    EXPECT_FALSE(bam.HasMergeQV());
    EXPECT_FALSE(bam.HasSubstitutionQV());
    EXPECT_FALSE(bam.HasSubstitutionTag());

    EXPECT_EQ(emptyString, bam.MovieName());
    EXPECT_EQ(emptyString, bam.ReadGroupId());

    // raw data
    tests::CheckRawData(bam);
}

TEST(BamRecordTest, FromBamRecordImpl)
{
    // check generic data
    BamRecordImpl genericBam = tests::CreateBamImpl();

    EXPECT_EQ(42, genericBam.Bin());
    EXPECT_EQ(42, genericBam.Flag());
    EXPECT_EQ(42, genericBam.InsertSize());
    EXPECT_EQ(42, genericBam.MapQuality());
    EXPECT_EQ(42, genericBam.MateReferenceId());
    EXPECT_EQ(42, genericBam.MatePosition());
    EXPECT_EQ(42, genericBam.Position());
    EXPECT_EQ(42, genericBam.ReferenceId());

    const TagCollection& genericTags = genericBam.Tags();
    EXPECT_TRUE(genericTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), genericTags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), genericTags.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), genericTags.at("CA").ToUInt8Array());

    // copy ctor
    BamRecord bam1(genericBam);

    EXPECT_EQ(42, bam1.impl_.Bin());
    EXPECT_EQ(42, bam1.impl_.Flag());
    EXPECT_EQ(42, bam1.impl_.InsertSize());
    EXPECT_EQ(42, bam1.impl_.MapQuality());
    EXPECT_EQ(42, bam1.impl_.MateReferenceId());
    EXPECT_EQ(42, bam1.impl_.MatePosition());
    EXPECT_EQ(42, bam1.impl_.Position());
    EXPECT_EQ(42, bam1.impl_.ReferenceId());

    const TagCollection& bam1Tags = bam1.impl_.Tags();
    EXPECT_TRUE(bam1Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), bam1Tags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), bam1Tags.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), bam1Tags.at("CA").ToUInt8Array());

    // copy assignment
    BamRecord bam2;
    bam2 = genericBam;

    EXPECT_EQ(42, bam2.impl_.Bin());
    EXPECT_EQ(42, bam2.impl_.Flag());
    EXPECT_EQ(42, bam2.impl_.InsertSize());
    EXPECT_EQ(42, bam2.impl_.MapQuality());
    EXPECT_EQ(42, bam2.impl_.MateReferenceId());
    EXPECT_EQ(42, bam2.impl_.MatePosition());
    EXPECT_EQ(42, bam2.impl_.Position());
    EXPECT_EQ(42, bam2.impl_.ReferenceId());

    const TagCollection& bam2Tags = bam2.impl_.Tags();
    EXPECT_TRUE(bam2Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), bam2Tags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), bam2Tags.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), bam2Tags.at("CA").ToUInt8Array());

    // change genericBam, make sure we deep copied bam1 & bam2
    genericBam.Position(2000);

    EXPECT_EQ(2000, genericBam.Position());
    EXPECT_EQ(42, bam1.impl_.Position());
    EXPECT_EQ(42, bam2.impl_.Position());

    // move ctor
    BamRecord bam3(move(tests::CreateBamImpl()));

    EXPECT_EQ(42, bam3.impl_.Bin());
    EXPECT_EQ(42, bam3.impl_.Flag());
    EXPECT_EQ(42, bam3.impl_.InsertSize());
    EXPECT_EQ(42, bam3.impl_.MapQuality());
    EXPECT_EQ(42, bam3.impl_.MateReferenceId());
    EXPECT_EQ(42, bam3.impl_.MatePosition());
    EXPECT_EQ(42, bam3.impl_.Position());
    EXPECT_EQ(42, bam3.impl_.ReferenceId());

    const TagCollection& bam3Tags = bam3.impl_.Tags();
    EXPECT_TRUE(bam3Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), bam3Tags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), bam3Tags.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), bam3Tags.at("CA").ToUInt8Array());

    // move assignment
    BamRecord bam4;
    bam4 = move(tests::CreateBamImpl());

    EXPECT_EQ(42, bam4.impl_.Bin());
    EXPECT_EQ(42, bam4.impl_.Flag());
    EXPECT_EQ(42, bam4.impl_.InsertSize());
    EXPECT_EQ(42, bam4.impl_.MapQuality());
    EXPECT_EQ(42, bam4.impl_.MateReferenceId());
    EXPECT_EQ(42, bam4.impl_.MatePosition());
    EXPECT_EQ(42, bam4.impl_.Position());
    EXPECT_EQ(42, bam4.impl_.ReferenceId());

    const TagCollection& bam4Tags = bam4.impl_.Tags();
    EXPECT_TRUE(bam4Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), bam4Tags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), bam4Tags.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), bam4Tags.at("CA").ToUInt8Array());
}

TEST(BamRecordTest, SelfAssignmentTolerated)
{
    BamRecord bam1;
    bam1.impl_.Bin(42);
    bam1.impl_.Flag(42);
    bam1.impl_.InsertSize(42);
    bam1.impl_.MapQuality(42);
    bam1.impl_.MatePosition(42);
    bam1.impl_.MateReferenceId(42);
    bam1.impl_.Position(42);
    bam1.impl_.ReferenceId(42);

    TagCollection tags;
    tags["HX"] = string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);
    bam1.impl_.Tags(tags);

    bam1 = bam1;

    EXPECT_EQ(42, bam1.impl_.Bin());
    EXPECT_EQ(42, bam1.impl_.Flag());
    EXPECT_EQ(42, bam1.impl_.InsertSize());
    EXPECT_EQ(42, bam1.impl_.MapQuality());
    EXPECT_EQ(42, bam1.impl_.MateReferenceId());
    EXPECT_EQ(42, bam1.impl_.MatePosition());
    EXPECT_EQ(42, bam1.impl_.Position());
    EXPECT_EQ(42, bam1.impl_.ReferenceId());

    const TagCollection& fetchedTags1 = bam1.impl_.Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    tests::CheckRawData(bam1);
}

TEST(BamRecordTest, CoreSetters)
{
    // create basic BAM with (generic) data
    BamRecord bam = tests::CreateBam();

    QualityValues testQVs;
    testQVs.push_back(0);
    testQVs.push_back(1);

    const string testTags = "GATTACA";

    // now set PacBio data
//    bam.AlignedStart(42);
//    bam.AlignedEnd(42);
//    bam.DeletionQVs(testQVs);
//    bam.DeletionTags(testTags);
//    bam.HoleNumber(42);
//    bam.InsertionQVs(testQVs);
//    bam.MergeQVs(testQVs);
//    bam.NumPasses(42);
//    bam.QueryEnd(42);
//    bam.QueryStart(42);
//    bam.ReadAccuracy(42);
//    bam.ReferenceEnd(42);
//    bam.ReferenceStart(42);
//    bam.SubstitutionQVs(testQVs);
//    bam.SubstitutionTags(testTags);

    // check generic data
    EXPECT_EQ(42, bam.impl_.Bin());
    EXPECT_EQ(42, bam.impl_.Flag());
    EXPECT_EQ(42, bam.impl_.InsertSize());
    EXPECT_EQ(42, bam.impl_.MapQuality());
    EXPECT_EQ(42, bam.impl_.MateReferenceId());
    EXPECT_EQ(42, bam.impl_.MatePosition());
    EXPECT_EQ(42, bam.impl_.Position());
    EXPECT_EQ(42, bam.impl_.ReferenceId());

    // check PacBio data
//    EXPECT_EQ(42, bam.AlignedStart());
//    EXPECT_EQ(42, bam.AlignedEnd());
//    EXPECT_EQ(testQVs, bam.DeletionQVs());
//    EXPECT_EQ(testTags, bam.DeletionTags());
//    EXPECT_EQ(42, bam.HoleNumber());
//    EXPECT_EQ(testQVs, bam.InsertionQVs());
//    EXPECT_EQ(testQVs, bam.MergeQVs());
//    EXPECT_EQ(42, bam.NumPasses());
//    EXPECT_EQ(42, bam.QueryEnd());
//    EXPECT_EQ(42, bam.QueryStart());
//    EXPECT_EQ(42, bam.ReadAccuracy());
//    EXPECT_EQ(42, bam.ReferenceEnd());
//    EXPECT_EQ(42, bam.ReferenceStart());
//    EXPECT_EQ(testQVs, bam.SubstitutionQVs());
//    EXPECT_EQ(testTags, bam.SubstitutionTags());

    // check tags
    const TagCollection& fetchedTags = bam.impl_.Tags();
    EXPECT_TRUE(fetchedTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(string("1abc75"), fetchedTags.at("HX").ToString());
    EXPECT_EQ(static_cast<int32_t>(-42), fetchedTags.at("XY").ToInt32());
    EXPECT_EQ(vector<uint8_t>({34, 5, 125}), fetchedTags.at("CA").ToUInt8Array());

    tests::CheckRawData(bam);
}

TEST(BamRecordTest, SequenceOrientation)
{
    const string sequence = "ATATATCCCGGCG";
    const string revSeq   = "CGCCGGGATATAT";

    // ----------------
    // forward strand
    // ----------------

    BamRecordImpl forwardImpl;
    forwardImpl.SetSequenceAndQualities(sequence);
    forwardImpl.SetReverseStrand(false);

    BamRecord forwardRead(forwardImpl);

    //  - "native" == "genomic"
    EXPECT_EQ(forwardRead.Sequence(Orientation::NATIVE),
              forwardRead.Sequence(Orientation::GENOMIC));
    //  - genomic output == genomic input
    EXPECT_EQ(sequence, forwardRead.Sequence(Orientation::GENOMIC));
    //  - native output == genomic input
    EXPECT_EQ(sequence, forwardRead.Sequence(Orientation::NATIVE));

    // ----------------
    // reverse strand
    // ----------------

    BamRecordImpl reverseImpl;
    reverseImpl.SetSequenceAndQualities(sequence);
    reverseImpl.SetReverseStrand(true);

    BamRecord reverseRead(reverseImpl);

    //  - "native" != "genomic"
    EXPECT_NE(reverseRead.Sequence(Orientation::NATIVE),
              reverseRead.Sequence(Orientation::GENOMIC));
    //  - genomic output == genomic input
    EXPECT_EQ(sequence, reverseRead.Sequence(Orientation::GENOMIC));
    //  - genomic (raw) input != native output
    EXPECT_NE(sequence, reverseRead.Sequence(Orientation::NATIVE));
    //  - native output should be reverse complement
    EXPECT_EQ(revSeq, reverseRead.Sequence(Orientation::NATIVE));
}

TEST(BamRecordTest, QualitiesOrientation)
{
    const string sequence = "ATATATCCCGGCG";
    const string qualities = "?]?]?]?]?]?]*";
    const string revQuals  = "*]?]?]?]?]?]?";

    // ----------------
    // forward strand
    // ----------------

    BamRecordImpl forwardImpl;
    forwardImpl.SetSequenceAndQualities(sequence, qualities);
    forwardImpl.SetReverseStrand(false);

    BamRecord forwardRead(forwardImpl);

    //  - "native" == "genomic"
    EXPECT_EQ(forwardRead.Qualities(Orientation::NATIVE),
              forwardRead.Qualities(Orientation::GENOMIC));
    //  - genomic (raw) input == genomic input
    EXPECT_EQ(qualities, forwardRead.Qualities(Orientation::GENOMIC).Fastq());
    //  - native output == genomic input
    EXPECT_EQ(qualities, forwardRead.Qualities(Orientation::NATIVE).Fastq());

    // ----------------
    // reverse strand
    // ----------------

    BamRecordImpl reverseImpl;
    reverseImpl.SetSequenceAndQualities(sequence, qualities);
    reverseImpl.SetReverseStrand(true);

    BamRecord reverseRead(reverseImpl);

    //  - "native" != "genomic"
    EXPECT_NE(reverseRead.Qualities(Orientation::NATIVE),
              reverseRead.Qualities(Orientation::GENOMIC));
    //  - genomic output == genomic input
    EXPECT_EQ(qualities, reverseRead.Qualities(Orientation::GENOMIC).Fastq());
    //  - genomic (raw) input != native output
    EXPECT_NE(qualities, reverseRead.Qualities(Orientation::NATIVE).Fastq());
    //  - native output should be reverse
    EXPECT_EQ(revQuals, reverseRead.Qualities(Orientation::NATIVE).Fastq());
}

TEST(BamRecordTest, SequenceTagsOrientation)
{
    const string tag    = "ATATATCCCGGCG";
    const string revTag = "CGCCGGGATATAT";

    // ----------------
    // forward strand
    // ----------------

    BamRecordImpl forwardImpl;
    forwardImpl.SetReverseStrand(false);
    forwardImpl.AddTag("dt", tag);
    forwardImpl.AddTag("st", tag);

    BamRecord forwardRead(forwardImpl);

    // sanity check
    EXPECT_TRUE(forwardImpl.HasTag("dt"));
    EXPECT_TRUE(forwardImpl.HasTag("st"));
    EXPECT_TRUE(forwardRead.HasDeletionTag());
    EXPECT_TRUE(forwardRead.HasSubstitutionTag());

    //  - "native" == "genomic"
    EXPECT_EQ(forwardRead.DeletionTag(Orientation::NATIVE),
              forwardRead.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(forwardRead.SubstitutionTag(Orientation::NATIVE),
              forwardRead.SubstitutionTag(Orientation::GENOMIC));

    //  - genomic output == genomic input
    EXPECT_EQ(tag, forwardRead.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tag, forwardRead.SubstitutionTag(Orientation::GENOMIC));

    //  - native output == genomic input
    EXPECT_EQ(tag, forwardRead.DeletionTag(Orientation::NATIVE));
    EXPECT_EQ(tag, forwardRead.SubstitutionTag(Orientation::NATIVE));

    // ----------------
    // reverse strand
    // ----------------

    BamRecordImpl reverseImpl;
    reverseImpl.SetReverseStrand(true);
    reverseImpl.AddTag("dt", revTag);
    reverseImpl.AddTag("st", revTag);

    BamRecord reverseRead(reverseImpl);

    // sanity check
    EXPECT_TRUE(reverseImpl.HasTag("dt"));
    EXPECT_TRUE(reverseImpl.HasTag("st"));
    EXPECT_TRUE(reverseRead.HasDeletionTag());
    EXPECT_TRUE(reverseRead.HasSubstitutionTag());

    //  - "native" != "genomic"
    EXPECT_NE(reverseRead.DeletionTag(Orientation::NATIVE),
              reverseRead.DeletionTag(Orientation::GENOMIC));
    EXPECT_NE(reverseRead.SubstitutionTag(Orientation::NATIVE),
              reverseRead.SubstitutionTag(Orientation::GENOMIC));

    //  - genomic output == genomic input
    EXPECT_EQ(tag, reverseRead.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(tag, reverseRead.SubstitutionTag(Orientation::GENOMIC));

    //  - genomic (raw) input != native output
    EXPECT_NE(tag, reverseRead.DeletionTag(Orientation::NATIVE));
    EXPECT_NE(tag, reverseRead.SubstitutionTag(Orientation::NATIVE));

    //  - native output should be reverse
    EXPECT_EQ(revTag, reverseRead.DeletionTag(Orientation::NATIVE));
    EXPECT_EQ(revTag, reverseRead.SubstitutionTag(Orientation::NATIVE));
}

TEST(BamRecordTest, FrameTagsOrientation)
{
    vector<uint16_t> frameData;
    vector<uint16_t> revFrameData;
    for (int i = 0, j = 4; i < 5; ++i, --j) {
        frameData.push_back(i*10);
        revFrameData.push_back(j*10);
    }
    const Frames frames(frameData);
    const Frames revFrames(revFrameData);

    // ----------------
    // forward strand
    // ----------------

    BamRecordImpl forwardImpl;
    forwardImpl.SetReverseStrand(false);
    forwardImpl.AddTag("ip", frames.Data());
    forwardImpl.AddTag("pw", frames.Data());

    BamRecord forwardRead(forwardImpl);

    // sanity check
    EXPECT_TRUE(forwardImpl.HasTag("ip"));
    EXPECT_TRUE(forwardImpl.HasTag("pw"));
    EXPECT_TRUE(forwardRead.HasIPD());
    EXPECT_TRUE(forwardRead.HasPulseWidth());

    //  - "native" == "genomic"
    EXPECT_EQ(forwardRead.IPD(Orientation::NATIVE),
              forwardRead.IPD(Orientation::GENOMIC));
    EXPECT_EQ(forwardRead.PulseWidth(Orientation::NATIVE),
              forwardRead.PulseWidth(Orientation::GENOMIC));

    //  - genomic output == genomic input
    EXPECT_EQ(frames, forwardRead.IPD(Orientation::GENOMIC));
    EXPECT_EQ(frames, forwardRead.PulseWidth(Orientation::GENOMIC));

    //  - native output == genomic input
    EXPECT_EQ(frames, forwardRead.IPD(Orientation::NATIVE));
    EXPECT_EQ(frames, forwardRead.PulseWidth(Orientation::NATIVE));

    // ----------------
    // reverse strand
    // ----------------

    BamRecordImpl reverseImpl;
    reverseImpl.SetReverseStrand(true);
    reverseImpl.AddTag("ip", revFrames.Data());
    reverseImpl.AddTag("pw", revFrames.Data());

    BamRecord reverseRead(reverseImpl);

    // sanity check
    EXPECT_TRUE(reverseImpl.HasTag("ip"));
    EXPECT_TRUE(reverseImpl.HasTag("pw"));
    EXPECT_TRUE(reverseRead.HasIPD());
    EXPECT_TRUE(reverseRead.HasPulseWidth());

    //  - "native" != "genomic"
    EXPECT_NE(reverseRead.IPD(Orientation::NATIVE),
              reverseRead.IPD(Orientation::GENOMIC));
    EXPECT_NE(reverseRead.PulseWidth(Orientation::NATIVE),
              reverseRead.PulseWidth(Orientation::GENOMIC));

    //  - genomic output == genomic input
    EXPECT_EQ(frames, reverseRead.IPD(Orientation::GENOMIC));
    EXPECT_EQ(frames, reverseRead.PulseWidth(Orientation::GENOMIC));

    //  - genomic (raw) input != native output
    EXPECT_NE(frames, reverseRead.IPD(Orientation::NATIVE));
    EXPECT_NE(frames, reverseRead.PulseWidth(Orientation::NATIVE));

    //  - native output should be reverse
    EXPECT_EQ(revFrames, reverseRead.IPD(Orientation::NATIVE));
    EXPECT_EQ(revFrames, reverseRead.PulseWidth(Orientation::NATIVE));
}

TEST(BamRecordTest, QualityTagsOrientation)
{
    const string qualities = "?]?]?]?]?]?]*";
    const string revQuals  = "*]?]?]?]?]?]?";

    // ----------------
    // forward strand
    // ----------------

    BamRecordImpl forwardImpl;
    forwardImpl.SetReverseStrand(false);
    forwardImpl.AddTag("dq", qualities);
    forwardImpl.AddTag("iq", qualities);
    forwardImpl.AddTag("mq", qualities);
    forwardImpl.AddTag("sq", qualities);

    BamRecord forwardRead(forwardImpl);

    // sanity check
    EXPECT_TRUE(forwardImpl.HasTag("dq"));
    EXPECT_TRUE(forwardImpl.HasTag("iq"));
    EXPECT_TRUE(forwardImpl.HasTag("mq"));
    EXPECT_TRUE(forwardImpl.HasTag("sq"));
    EXPECT_TRUE(forwardRead.HasDeletionQV());
    EXPECT_TRUE(forwardRead.HasInsertionQV());
    EXPECT_TRUE(forwardRead.HasMergeQV());
    EXPECT_TRUE(forwardRead.HasSubstitutionQV());

    //  - "native" == "genomic"
    EXPECT_EQ(forwardRead.DeletionQV(Orientation::NATIVE),
              forwardRead.DeletionQV(Orientation::GENOMIC));
    EXPECT_EQ(forwardRead.InsertionQV(Orientation::NATIVE),
              forwardRead.InsertionQV(Orientation::GENOMIC));
    EXPECT_EQ(forwardRead.MergeQV(Orientation::NATIVE),
              forwardRead.MergeQV(Orientation::GENOMIC));
    EXPECT_EQ(forwardRead.SubstitutionQV(Orientation::NATIVE),
              forwardRead.SubstitutionQV(Orientation::GENOMIC));

    //  - genomic output == genomic input
    EXPECT_EQ(qualities, forwardRead.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(qualities, forwardRead.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(qualities, forwardRead.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(qualities, forwardRead.SubstitutionQV(Orientation::GENOMIC).Fastq());

    //  - native output == genomic input
    EXPECT_EQ(qualities, forwardRead.DeletionQV(Orientation::NATIVE).Fastq());
    EXPECT_EQ(qualities, forwardRead.InsertionQV(Orientation::NATIVE).Fastq());
    EXPECT_EQ(qualities, forwardRead.MergeQV(Orientation::NATIVE).Fastq());
    EXPECT_EQ(qualities, forwardRead.SubstitutionQV(Orientation::NATIVE).Fastq());

    // ----------------
    // reverse strand
    // ----------------

    BamRecordImpl reverseImpl;
    reverseImpl.SetReverseStrand(true);
    reverseImpl.AddTag("dq", revQuals);
    reverseImpl.AddTag("iq", revQuals);
    reverseImpl.AddTag("mq", revQuals);
    reverseImpl.AddTag("sq", revQuals);

    BamRecord reverseRead(reverseImpl);

    // sanity check
    EXPECT_TRUE(reverseImpl.HasTag("dq"));
    EXPECT_TRUE(reverseImpl.HasTag("iq"));
    EXPECT_TRUE(reverseImpl.HasTag("mq"));
    EXPECT_TRUE(reverseImpl.HasTag("sq"));
    EXPECT_TRUE(reverseRead.HasDeletionQV());
    EXPECT_TRUE(reverseRead.HasInsertionQV());
    EXPECT_TRUE(reverseRead.HasMergeQV());
    EXPECT_TRUE(reverseRead.HasSubstitutionQV());

    //  - "native" != "genomic"
    EXPECT_NE(reverseRead.DeletionQV(Orientation::NATIVE),
              reverseRead.DeletionQV(Orientation::GENOMIC));
    EXPECT_NE(reverseRead.InsertionQV(Orientation::NATIVE),
              reverseRead.InsertionQV(Orientation::GENOMIC));
    EXPECT_NE(reverseRead.MergeQV(Orientation::NATIVE),
              reverseRead.MergeQV(Orientation::GENOMIC));
    EXPECT_NE(reverseRead.SubstitutionQV(Orientation::NATIVE),
              reverseRead.SubstitutionQV(Orientation::GENOMIC));

    //  - genomic output == genomic input
    EXPECT_EQ(qualities, reverseRead.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(qualities, reverseRead.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(qualities, reverseRead.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(qualities, reverseRead.SubstitutionQV(Orientation::GENOMIC).Fastq());

    //  - genomic (raw) input != native output
    EXPECT_NE(qualities, reverseRead.DeletionQV(Orientation::NATIVE).Fastq());
    EXPECT_NE(qualities, reverseRead.InsertionQV(Orientation::NATIVE).Fastq());
    EXPECT_NE(qualities, reverseRead.MergeQV(Orientation::NATIVE).Fastq());
    EXPECT_NE(qualities, reverseRead.SubstitutionQV(Orientation::NATIVE).Fastq());

    //  - native output should be reverse
    EXPECT_EQ(revQuals, reverseRead.DeletionQV(Orientation::NATIVE).Fastq());
    EXPECT_EQ(revQuals, reverseRead.InsertionQV(Orientation::NATIVE).Fastq());
    EXPECT_EQ(revQuals, reverseRead.MergeQV(Orientation::NATIVE).Fastq());
    EXPECT_EQ(revQuals, reverseRead.SubstitutionQV(Orientation::NATIVE).Fastq());
}

TEST(BamRecordTest, ClippingAndOrientation)
{
    // forward string names, cigar, seq
    // reverse strand records have same cigar and **input** seq as forward strand
    // (native output will be rev-comp'd)

    const string s1_cigar = "10M";
    const string s2_cigar = "3M4N3M";
    const string s3_cigar = "1S8M1S";
    const string s4_cigar = "1H8M1H";
    const string s5_cigar = "2S6M2S";
    const string s6_cigar = "2S3M2I3M2S";
    const string s7_cigar = "2H6M2H";

    const string s1_seq  = "ATCCGCGGTT";
    const string s2_seq  = "ACGTT";
    const string s3_seq  = "ACCCGCGGTT";
    const string s4_seq  = "ATCGCGGT";
    const string s5_seq  = "AGCCGCGGTT";
    const string s6_seq  = "ATCCGNNCGGTT";
    const string s7_seq  = "CAGCGG";

    const string s1_seq_clipped  = "ATCCGCGGTT";
    const string s2_seq_clipped  = "ACGTT";
    const string s3_seq_clipped  = "CCCGCGGT";
    const string s4_seq_clipped  = "ATCGCGGT";
    const string s5_seq_clipped  = "CCGCGG";
    const string s6_seq_clipped  = "CCGNNCGG";
    const string s7_seq_clipped  = "CAGCGG";

    const string s1_revseq = "AACCGCGGAT";
    const string s2_revseq = "AACGT";
    const string s3_revseq = "AACCGCGGGT";
    const string s4_revseq = "ACCGCGAT";
    const string s5_revseq = "AACCGCGGCT";
    const string s6_revseq = "AACCGNNCGGAT";
    const string s7_revseq = "CCGCTG";

    const string s1_revseq_clipped = "AACCGCGGAT";
    const string s2_revseq_clipped = "AACGT";
    const string s3_revseq_clipped = "ACCGCGGG";
    const string s4_revseq_clipped = "ACCGCGAT";
    const string s5_revseq_clipped = "CCGCGG";
    const string s6_revseq_clipped = "CCGNNCGG";
    const string s7_revseq_clipped = "CCGCTG";

    const BamRecord s1 = tests::MakeCigaredRecord(s1_seq, s1_cigar, false);
    const BamRecord s2 = tests::MakeCigaredRecord(s2_seq, s2_cigar, false);
    const BamRecord s3 = tests::MakeCigaredRecord(s3_seq, s3_cigar, false);
    const BamRecord s4 = tests::MakeCigaredRecord(s4_seq, s4_cigar, false);
    const BamRecord s5 = tests::MakeCigaredRecord(s5_seq, s5_cigar, false);
    const BamRecord s6 = tests::MakeCigaredRecord(s6_seq, s6_cigar, false);
    const BamRecord s7 = tests::MakeCigaredRecord(s7_seq, s7_cigar, false);
    const BamRecord s1_reverse = tests::MakeCigaredRecord(s1_seq, s1_cigar, true);
    const BamRecord s2_reverse = tests::MakeCigaredRecord(s2_seq, s2_cigar, true);
    const BamRecord s3_reverse = tests::MakeCigaredRecord(s3_seq, s3_cigar, true);
    const BamRecord s4_reverse = tests::MakeCigaredRecord(s4_seq, s4_cigar, true);
    const BamRecord s5_reverse = tests::MakeCigaredRecord(s5_seq, s5_cigar, true);
    const BamRecord s6_reverse = tests::MakeCigaredRecord(s6_seq, s6_cigar, true);
    const BamRecord s7_reverse = tests::MakeCigaredRecord(s7_seq, s7_cigar, true);

    // ----------------
    // forward strand
    // ----------------

    //  - "native" == "genomic"
    EXPECT_EQ(s1.Sequence(Orientation::NATIVE), s1.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2.Sequence(Orientation::NATIVE), s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3.Sequence(Orientation::NATIVE), s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s4.Sequence(Orientation::NATIVE), s4.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s5.Sequence(Orientation::NATIVE), s5.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s6.Sequence(Orientation::NATIVE), s6.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s7.Sequence(Orientation::NATIVE), s7.Sequence(Orientation::GENOMIC));

    //  - unclipped, unaligned genomic output == genomic input
    EXPECT_EQ(s1_seq, s1.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_seq, s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_seq, s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s4_seq, s4.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s5_seq, s5.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s6_seq, s6.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s7_seq, s7.Sequence(Orientation::GENOMIC));

    //  - unclipped, unaligned native output == genomic input
    EXPECT_EQ(s1_seq, s1.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s2_seq, s2.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s3_seq, s3.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s4_seq, s4.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s5_seq, s5.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s6_seq, s6.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s7_seq, s7.Sequence(Orientation::NATIVE));

    //  - clipped, unaligned native output == clipped genomic input
    EXPECT_EQ(s1_seq_clipped, s1.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s2_seq_clipped, s2.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s3_seq_clipped, s3.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s4_seq_clipped, s4.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s5_seq_clipped, s5.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s6_seq_clipped, s6.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s7_seq_clipped, s7.Sequence(Orientation::NATIVE, false, true));

    // ----------------
    // reverse strand
    // ----------------

    //  - "native" != "genomic"
    EXPECT_NE(s1_reverse.Sequence(Orientation::NATIVE), s1_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.Sequence(Orientation::NATIVE), s2_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.Sequence(Orientation::NATIVE), s3_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.Sequence(Orientation::NATIVE), s4_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.Sequence(Orientation::NATIVE), s5_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.Sequence(Orientation::NATIVE), s6_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s7_reverse.Sequence(Orientation::NATIVE), s7_reverse.Sequence(Orientation::GENOMIC));

    //  - unclipped, unaligned genomic output == genomic input
    EXPECT_EQ(s1_seq, s1_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2_seq, s2_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3_seq, s3_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s4_seq, s4_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s5_seq, s5_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s6_seq, s6_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s7_seq, s7_reverse.Sequence(Orientation::GENOMIC));

    //  - unclipped, unaligned native output != genomic (raw) input
    EXPECT_NE(s1_seq, s1_reverse.Sequence(Orientation::NATIVE));
    EXPECT_NE(s2_seq, s2_reverse.Sequence(Orientation::NATIVE));
    EXPECT_NE(s3_seq, s3_reverse.Sequence(Orientation::NATIVE));
    EXPECT_NE(s4_seq, s4_reverse.Sequence(Orientation::NATIVE));
    EXPECT_NE(s5_seq, s5_reverse.Sequence(Orientation::NATIVE));
    EXPECT_NE(s6_seq, s6_reverse.Sequence(Orientation::NATIVE));
    EXPECT_NE(s7_seq, s7_reverse.Sequence(Orientation::NATIVE));

    //  - unclipped, unaligned native output should be reverse
    EXPECT_EQ(s1_revseq, s1_reverse.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s2_revseq, s2_reverse.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s3_revseq, s3_reverse.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s4_revseq, s4_reverse.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s5_revseq, s5_reverse.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s6_revseq, s6_reverse.Sequence(Orientation::NATIVE));
    EXPECT_EQ(s7_revseq, s7_reverse.Sequence(Orientation::NATIVE));

    //  - clipped, unaligned native output == clipped genomic input
    EXPECT_EQ(s1_revseq_clipped, s1_reverse.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s2_revseq_clipped, s2_reverse.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s3_revseq_clipped, s3_reverse.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s4_revseq_clipped, s4_reverse.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s5_revseq_clipped, s5_reverse.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s6_revseq_clipped, s6_reverse.Sequence(Orientation::NATIVE, false, true));
    EXPECT_EQ(s7_revseq_clipped, s7_reverse.Sequence(Orientation::NATIVE, false, true));
}

TEST(BamRecordTest, ClippingOrientationAndAlignment)
{
    // forward string names, cigar, seq
    // reverse strand records have same cigar and **input** seq as forward strand
    // (native output will be rev-comp'd)

    const string s1_cigar = "4M3D4M";
    const string s2_cigar = "4M1D2I2D4M";
    const string s3_cigar = "4M1D2P2I2P2D4M";
    const string s4_cigar = "2S4M3D4M3S";
    const string s5_cigar = "2H4M3D4M3H";
    const string s6_cigar = "2H2S4M3D4M3S3H";

    const string s1_seq = "AACCGTTA";
    const string s2_seq = "ATCCTAGGTT";
    const string s3_seq = "ATCCTAGGTT";
    const string s4_seq = "TTAACCGTTACCG";
    const string s5_seq = "AACCGTTA";
    const string s6_seq = "TTAACCGTTACCG";

    const string s1_seq_aligned = "AACC---GTTA";
    const string s2_seq_aligned = "ATCC-TA--GGTT";
    const string s3_seq_aligned = "ATCC-**TA**--GGTT";
    const string s4_seq_aligned = "TTAACC---GTTACCG";
    const string s5_seq_aligned = "AACC---GTTA";
    const string s6_seq_aligned = "TTAACC---GTTACCG";

    const string s1_seq_aligned_clipped = "AACC---GTTA";
    const string s2_seq_aligned_clipped = "ATCC-TA--GGTT";
    const string s3_seq_aligned_clipped = "ATCC-**TA**--GGTT";
    const string s4_seq_aligned_clipped = "AACC---GTTA";
    const string s5_seq_aligned_clipped = "AACC---GTTA";
    const string s6_seq_aligned_clipped = "AACC---GTTA";

    const string s1_revseq = "TAACGGTT";
    const string s2_revseq = "AACCTAGGAT";
    const string s3_revseq = "AACCTAGGAT";
    const string s4_revseq = "CGGTAACGGTTAA";
    const string s5_revseq = "TAACGGTT";
    const string s6_revseq = "CGGTAACGGTTAA";

    const string s1_revseq_aligned = "TAAC---GGTT";
    const string s2_revseq_aligned = "AACC--TA-GGAT";
    const string s3_revseq_aligned = "AACC--**TA**-GGAT";
    const string s4_revseq_aligned = "CGGTAAC---GGTTAA";
    const string s5_revseq_aligned = "TAAC---GGTT";
    const string s6_revseq_aligned = "CGGTAAC---GGTTAA";

    const string s1_revseq_aligned_clipped = "TAAC---GGTT";
    const string s2_revseq_aligned_clipped = "AACC--TA-GGAT";
    const string s3_revseq_aligned_clipped = "AACC--**TA**-GGAT";
    const string s4_revseq_aligned_clipped = "TAAC---GGTT";
    const string s5_revseq_aligned_clipped = "TAAC---GGTT";
    const string s6_revseq_aligned_clipped = "TAAC---GGTT";

    const BamRecord s1 = tests::MakeCigaredRecord(s1_seq, s1_cigar, false);
    const BamRecord s2 = tests::MakeCigaredRecord(s2_seq, s2_cigar, false);
    const BamRecord s3 = tests::MakeCigaredRecord(s3_seq, s3_cigar, false);
    const BamRecord s4 = tests::MakeCigaredRecord(s4_seq, s4_cigar, false);
    const BamRecord s5 = tests::MakeCigaredRecord(s5_seq, s5_cigar, false);
    const BamRecord s6 = tests::MakeCigaredRecord(s6_seq, s6_cigar, false);
    const BamRecord s1_reverse = tests::MakeCigaredRecord(s1_seq, s1_cigar, true);
    const BamRecord s2_reverse = tests::MakeCigaredRecord(s2_seq, s2_cigar, true);
    const BamRecord s3_reverse = tests::MakeCigaredRecord(s3_seq, s3_cigar, true);
    const BamRecord s4_reverse = tests::MakeCigaredRecord(s4_seq, s4_cigar, true);
    const BamRecord s5_reverse = tests::MakeCigaredRecord(s5_seq, s5_cigar, true);
    const BamRecord s6_reverse = tests::MakeCigaredRecord(s6_seq, s6_cigar, true);

    // ----------------
    // forward strand
    // ----------------

    //  - "native" == "genomic"
    EXPECT_EQ(s1.Sequence(Orientation::NATIVE), s1.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s2.Sequence(Orientation::NATIVE), s2.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s3.Sequence(Orientation::NATIVE), s3.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s4.Sequence(Orientation::NATIVE), s4.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s5.Sequence(Orientation::NATIVE), s5.Sequence(Orientation::GENOMIC));
    EXPECT_EQ(s6.Sequence(Orientation::NATIVE), s6.Sequence(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output == genomic input
    EXPECT_EQ(s1_seq, s1.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s2_seq, s2.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s3_seq, s3.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s4_seq, s4.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s5_seq, s5.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s6_seq, s6.Sequence(Orientation::GENOMIC, false, false));

    //  - aligned, unclipped genomic output == aligned, unclipped genomic input
    EXPECT_EQ(s1_seq_aligned, s1.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s2_seq_aligned, s2.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s3_seq_aligned, s3.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s4_seq_aligned, s4.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s5_seq_aligned, s5.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s6_seq_aligned, s6.Sequence(Orientation::GENOMIC, true, false));

    //  - aligned, clipped genomic output == aligned, clipped genomic output
    EXPECT_EQ(s1_seq_aligned_clipped, s1.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s2_seq_aligned_clipped, s2.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s3_seq_aligned_clipped, s3.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s4_seq_aligned_clipped, s4.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s5_seq_aligned_clipped, s5.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s6_seq_aligned_clipped, s6.Sequence(Orientation::GENOMIC, true, true));

    // ----------------
    // reverse strand
    // ----------------

    //  - "native" != "genomic"
    EXPECT_NE(s1_reverse.Sequence(Orientation::NATIVE), s1_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.Sequence(Orientation::NATIVE), s2_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.Sequence(Orientation::NATIVE), s3_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.Sequence(Orientation::NATIVE), s4_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.Sequence(Orientation::NATIVE), s5_reverse.Sequence(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.Sequence(Orientation::NATIVE), s6_reverse.Sequence(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output
    EXPECT_EQ(s1_seq, s1_reverse.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s2_seq, s2_reverse.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s3_seq, s3_reverse.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s4_seq, s4_reverse.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s5_seq, s5_reverse.Sequence(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s6_seq, s6_reverse.Sequence(Orientation::GENOMIC, false, false));

    //  - unaligned, unclipped native output
    EXPECT_EQ(s1_revseq, s1_reverse.Sequence(Orientation::NATIVE, false, false));
    EXPECT_EQ(s2_revseq, s2_reverse.Sequence(Orientation::NATIVE, false, false));
    EXPECT_EQ(s3_revseq, s3_reverse.Sequence(Orientation::NATIVE, false, false));
    EXPECT_EQ(s4_revseq, s4_reverse.Sequence(Orientation::NATIVE, false, false));
    EXPECT_EQ(s5_revseq, s5_reverse.Sequence(Orientation::NATIVE, false, false));
    EXPECT_EQ(s6_revseq, s6_reverse.Sequence(Orientation::NATIVE, false, false));

    //  - aligned, unclipped genomic output
    EXPECT_EQ(s1_seq_aligned, s1_reverse.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s2_seq_aligned, s2_reverse.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s3_seq_aligned, s3_reverse.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s4_seq_aligned, s4_reverse.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s5_seq_aligned, s5_reverse.Sequence(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s6_seq_aligned, s6_reverse.Sequence(Orientation::GENOMIC, true, false));

    //  - aligned, unclipped native output
    EXPECT_EQ(s1_revseq_aligned, s1_reverse.Sequence(Orientation::NATIVE, true, false));
    EXPECT_EQ(s2_revseq_aligned, s2_reverse.Sequence(Orientation::NATIVE, true, false));
    EXPECT_EQ(s3_revseq_aligned, s3_reverse.Sequence(Orientation::NATIVE, true, false));
    EXPECT_EQ(s4_revseq_aligned, s4_reverse.Sequence(Orientation::NATIVE, true, false));
    EXPECT_EQ(s5_revseq_aligned, s5_reverse.Sequence(Orientation::NATIVE, true, false));
    EXPECT_EQ(s6_revseq_aligned, s6_reverse.Sequence(Orientation::NATIVE, true, false));

    //  - aligned, clipped genomic output
    EXPECT_EQ(s1_seq_aligned_clipped, s1_reverse.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s2_seq_aligned_clipped, s2_reverse.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s3_seq_aligned_clipped, s3_reverse.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s4_seq_aligned_clipped, s4_reverse.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s5_seq_aligned_clipped, s5_reverse.Sequence(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s6_seq_aligned_clipped, s6_reverse.Sequence(Orientation::GENOMIC, true, true));

    //  - aligned, clipped native output
    EXPECT_EQ(s1_revseq_aligned_clipped, s1_reverse.Sequence(Orientation::NATIVE, true, true));
    EXPECT_EQ(s2_revseq_aligned_clipped, s2_reverse.Sequence(Orientation::NATIVE, true, true));
    EXPECT_EQ(s3_revseq_aligned_clipped, s3_reverse.Sequence(Orientation::NATIVE, true, true));
    EXPECT_EQ(s4_revseq_aligned_clipped, s4_reverse.Sequence(Orientation::NATIVE, true, true));
    EXPECT_EQ(s5_revseq_aligned_clipped, s5_reverse.Sequence(Orientation::NATIVE, true, true));
    EXPECT_EQ(s6_revseq_aligned_clipped, s6_reverse.Sequence(Orientation::NATIVE, true, true));
}


TEST(BamRecordTest, QualityTagsClippedAndAligned)
{
    // NOTE - FASTQ for QV=0 is '!'. Thus deletions/padding will be rendered as '!'s.

    const string s1_cigar = "4M3D4M";
    const string s2_cigar = "4M1D2I2D4M";
    const string s3_cigar = "4M1D2P2I2P2D4M";
    const string s4_cigar = "3S4M3D4M3S";
    const string s5_cigar = "2H4M3D4M3H";
    const string s6_cigar = "2H3S4M3D4M3S3H";

    const string s1_quals = "?]?]?]?@";
    const string s2_quals = "?]?]87?]?@";
    const string s3_quals = "?]?]87?]?@";
    const string s4_quals = "vvv?]?]?]?@xxx";
    const string s5_quals = "?]?]?]?@";
    const string s6_quals = "vvv?]?]?]?@xxx";

    const string s1_quals_aligned = "?]?]!!!?]?@";
    const string s2_quals_aligned = "?]?]!87!!?]?@";
    const string s3_quals_aligned = "?]?]!!!87!!!!?]?@";
    const string s4_quals_aligned = "vvv?]?]!!!?]?@xxx";
    const string s5_quals_aligned = "?]?]!!!?]?@";
    const string s6_quals_aligned = "vvv?]?]!!!?]?@xxx";

    const string s1_quals_aligned_clipped = "?]?]!!!?]?@";
    const string s2_quals_aligned_clipped = "?]?]!87!!?]?@";
    const string s3_quals_aligned_clipped = "?]?]!!!87!!!!?]?@";
    const string s4_quals_aligned_clipped = "?]?]!!!?]?@";
    const string s5_quals_aligned_clipped = "?]?]!!!?]?@";
    const string s6_quals_aligned_clipped = "?]?]!!!?]?@";

    const string s1_revquals = "@?]?]?]?";
    const string s2_revquals = "@?]?78]?]?";
    const string s3_revquals = "@?]?78]?]?";
    const string s4_revquals = "xxx@?]?]?]?vvv";
    const string s5_revquals = "@?]?]?]?";
    const string s6_revquals = "xxx@?]?]?]?vvv";

    const string s1_revquals_aligned = "@?]?!!!]?]?";
    const string s2_revquals_aligned = "@?]?!78!!]?]?";
    const string s3_revquals_aligned = "@?]?!!!78!!!!]?]?";
    const string s4_revquals_aligned = "xxx@?]?!!!]?]?vvv";
    const string s5_revquals_aligned = "@?]?!!!]?]?";
    const string s6_revquals_aligned = "xxx@?]?!!!]?]?vvv";

    const string s1_revquals_aligned_clipped = "@?]?!!!]?]?";
    const string s2_revquals_aligned_clipped = "@?]?!78!!]?]?";
    const string s3_revquals_aligned_clipped = "@?]?!!!78!!!!]?]?";
    const string s4_revquals_aligned_clipped = "@?]?!!!]?]?";
    const string s5_revquals_aligned_clipped = "@?]?!!!]?]?";
    const string s6_revquals_aligned_clipped = "@?]?!!!]?]?";

    const BamRecord s1 = tests::MakeCigaredQualRecord(s1_quals, s1_cigar, false);
    const BamRecord s2 = tests::MakeCigaredQualRecord(s2_quals, s2_cigar, false);
    const BamRecord s3 = tests::MakeCigaredQualRecord(s3_quals, s3_cigar, false);
    const BamRecord s4 = tests::MakeCigaredQualRecord(s4_quals, s4_cigar, false);
    const BamRecord s5 = tests::MakeCigaredQualRecord(s5_quals, s5_cigar, false);
    const BamRecord s6 = tests::MakeCigaredQualRecord(s6_quals, s6_cigar, false);
    const BamRecord s1_reverse = tests::MakeCigaredQualRecord(s1_quals, s1_cigar, true);
    const BamRecord s2_reverse = tests::MakeCigaredQualRecord(s2_quals, s2_cigar, true);
    const BamRecord s3_reverse = tests::MakeCigaredQualRecord(s3_quals, s3_cigar, true);
    const BamRecord s4_reverse = tests::MakeCigaredQualRecord(s4_quals, s4_cigar, true);
    const BamRecord s5_reverse = tests::MakeCigaredQualRecord(s5_quals, s5_cigar, true);
    const BamRecord s6_reverse = tests::MakeCigaredQualRecord(s6_quals, s6_cigar, true);

    // ----------------
    // forward strand
    // ----------------

    //  - "native" == "genomic"
    EXPECT_EQ(s1.DeletionQV(Orientation::NATIVE), s1.DeletionQV(Orientation::GENOMIC));
    EXPECT_EQ(s2.DeletionQV(Orientation::NATIVE), s2.DeletionQV(Orientation::GENOMIC));
    EXPECT_EQ(s3.DeletionQV(Orientation::NATIVE), s3.DeletionQV(Orientation::GENOMIC));
    EXPECT_EQ(s4.DeletionQV(Orientation::NATIVE), s4.DeletionQV(Orientation::GENOMIC));
    EXPECT_EQ(s5.DeletionQV(Orientation::NATIVE), s5.DeletionQV(Orientation::GENOMIC));
    EXPECT_EQ(s6.DeletionQV(Orientation::NATIVE), s6.DeletionQV(Orientation::GENOMIC));

    EXPECT_EQ(s1.InsertionQV(Orientation::NATIVE), s1.InsertionQV(Orientation::GENOMIC));
    EXPECT_EQ(s2.InsertionQV(Orientation::NATIVE), s2.InsertionQV(Orientation::GENOMIC));
    EXPECT_EQ(s3.InsertionQV(Orientation::NATIVE), s3.InsertionQV(Orientation::GENOMIC));
    EXPECT_EQ(s4.InsertionQV(Orientation::NATIVE), s4.InsertionQV(Orientation::GENOMIC));
    EXPECT_EQ(s5.InsertionQV(Orientation::NATIVE), s5.InsertionQV(Orientation::GENOMIC));
    EXPECT_EQ(s6.InsertionQV(Orientation::NATIVE), s6.InsertionQV(Orientation::GENOMIC));

    EXPECT_EQ(s1.MergeQV(Orientation::NATIVE), s1.MergeQV(Orientation::GENOMIC));
    EXPECT_EQ(s2.MergeQV(Orientation::NATIVE), s2.MergeQV(Orientation::GENOMIC));
    EXPECT_EQ(s3.MergeQV(Orientation::NATIVE), s3.MergeQV(Orientation::GENOMIC));
    EXPECT_EQ(s4.MergeQV(Orientation::NATIVE), s4.MergeQV(Orientation::GENOMIC));
    EXPECT_EQ(s5.MergeQV(Orientation::NATIVE), s5.MergeQV(Orientation::GENOMIC));
    EXPECT_EQ(s6.MergeQV(Orientation::NATIVE), s6.MergeQV(Orientation::GENOMIC));

    EXPECT_EQ(s1.SubstitutionQV(Orientation::NATIVE), s1.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_EQ(s2.SubstitutionQV(Orientation::NATIVE), s2.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_EQ(s3.SubstitutionQV(Orientation::NATIVE), s3.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_EQ(s4.SubstitutionQV(Orientation::NATIVE), s4.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_EQ(s5.SubstitutionQV(Orientation::NATIVE), s5.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_EQ(s6.SubstitutionQV(Orientation::NATIVE), s6.SubstitutionQV(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output == genomic input
    EXPECT_EQ(s1_quals, s1.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_quals, s2.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_quals, s3.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s4_quals, s4.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s5_quals, s5.DeletionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s6_quals, s6.DeletionQV(Orientation::GENOMIC).Fastq());

    EXPECT_EQ(s1_quals, s1.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_quals, s2.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_quals, s3.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s4_quals, s4.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s5_quals, s5.InsertionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s6_quals, s6.InsertionQV(Orientation::GENOMIC).Fastq());

    EXPECT_EQ(s1_quals, s1.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_quals, s2.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_quals, s3.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s4_quals, s4.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s5_quals, s5.MergeQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s6_quals, s6.MergeQV(Orientation::GENOMIC).Fastq());

    EXPECT_EQ(s1_quals, s1.SubstitutionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s2_quals, s2.SubstitutionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s3_quals, s3.SubstitutionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s4_quals, s4.SubstitutionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s5_quals, s5.SubstitutionQV(Orientation::GENOMIC).Fastq());
    EXPECT_EQ(s6_quals, s6.SubstitutionQV(Orientation::GENOMIC).Fastq());

    //  - aligned, unclipped genomic output == aligned, unclipped genomic input
    EXPECT_EQ(s1_quals_aligned, s1.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6.DeletionQV(Orientation::GENOMIC, true, false).Fastq());

    EXPECT_EQ(s1_quals_aligned, s1.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6.MergeQV(Orientation::GENOMIC, true, false).Fastq());

    EXPECT_EQ(s1_quals_aligned, s1.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6.InsertionQV(Orientation::GENOMIC, true, false).Fastq());

    EXPECT_EQ(s1_quals_aligned, s1.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());

    //  - aligned, clipped genomic output == aligned, clipped genomic output
    EXPECT_EQ(s1_quals_aligned_clipped, s1.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6.DeletionQV(Orientation::GENOMIC, true, true).Fastq());

    EXPECT_EQ(s1_quals_aligned_clipped, s1.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6.MergeQV(Orientation::GENOMIC, true, true).Fastq());

    EXPECT_EQ(s1_quals_aligned_clipped, s1.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6.InsertionQV(Orientation::GENOMIC, true, true).Fastq());

    EXPECT_EQ(s1_quals_aligned_clipped, s1.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());

    // ----------------
    // reverse strand
    // ----------------

    //  - "native" != "genomic"
    EXPECT_NE(s1_reverse.DeletionQV(Orientation::NATIVE), s1_reverse.DeletionQV(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.DeletionQV(Orientation::NATIVE), s2_reverse.DeletionQV(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.DeletionQV(Orientation::NATIVE), s3_reverse.DeletionQV(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.DeletionQV(Orientation::NATIVE), s4_reverse.DeletionQV(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.DeletionQV(Orientation::NATIVE), s5_reverse.DeletionQV(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.DeletionQV(Orientation::NATIVE), s6_reverse.DeletionQV(Orientation::GENOMIC));

    EXPECT_NE(s1_reverse.InsertionQV(Orientation::NATIVE), s1_reverse.InsertionQV(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.InsertionQV(Orientation::NATIVE), s2_reverse.InsertionQV(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.InsertionQV(Orientation::NATIVE), s3_reverse.InsertionQV(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.InsertionQV(Orientation::NATIVE), s4_reverse.InsertionQV(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.InsertionQV(Orientation::NATIVE), s5_reverse.InsertionQV(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.InsertionQV(Orientation::NATIVE), s6_reverse.InsertionQV(Orientation::GENOMIC));

    EXPECT_NE(s1_reverse.MergeQV(Orientation::NATIVE), s1_reverse.MergeQV(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.MergeQV(Orientation::NATIVE), s2_reverse.MergeQV(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.MergeQV(Orientation::NATIVE), s3_reverse.MergeQV(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.MergeQV(Orientation::NATIVE), s4_reverse.MergeQV(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.MergeQV(Orientation::NATIVE), s5_reverse.MergeQV(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.MergeQV(Orientation::NATIVE), s6_reverse.MergeQV(Orientation::GENOMIC));

    EXPECT_NE(s1_reverse.SubstitutionQV(Orientation::NATIVE), s1_reverse.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.SubstitutionQV(Orientation::NATIVE), s2_reverse.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.SubstitutionQV(Orientation::NATIVE), s3_reverse.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.SubstitutionQV(Orientation::NATIVE), s4_reverse.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.SubstitutionQV(Orientation::NATIVE), s5_reverse.SubstitutionQV(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.SubstitutionQV(Orientation::NATIVE), s6_reverse.SubstitutionQV(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revquals, s1_reverse.DeletionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s2_revquals, s2_reverse.DeletionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s3_revquals, s3_reverse.DeletionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s4_revquals, s4_reverse.DeletionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s5_revquals, s5_reverse.DeletionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s6_revquals, s6_reverse.DeletionQV(Orientation::GENOMIC, false, false).Fastq());

    EXPECT_EQ(s1_revquals, s1_reverse.InsertionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s2_revquals, s2_reverse.InsertionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s3_revquals, s3_reverse.InsertionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s4_revquals, s4_reverse.InsertionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s5_revquals, s5_reverse.InsertionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s6_revquals, s6_reverse.InsertionQV(Orientation::GENOMIC, false, false).Fastq());

    EXPECT_EQ(s1_revquals, s1_reverse.MergeQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s2_revquals, s2_reverse.MergeQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s3_revquals, s3_reverse.MergeQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s4_revquals, s4_reverse.MergeQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s5_revquals, s5_reverse.MergeQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s6_revquals, s6_reverse.MergeQV(Orientation::GENOMIC, false, false).Fastq());

    EXPECT_EQ(s1_revquals, s1_reverse.SubstitutionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s2_revquals, s2_reverse.SubstitutionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s3_revquals, s3_reverse.SubstitutionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s4_revquals, s4_reverse.SubstitutionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s5_revquals, s5_reverse.SubstitutionQV(Orientation::GENOMIC, false, false).Fastq());
    EXPECT_EQ(s6_revquals, s6_reverse.SubstitutionQV(Orientation::GENOMIC, false, false).Fastq());

    //  - unaligned, unclipped native output (native input)
    EXPECT_EQ(s1_quals, s1_reverse.DeletionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s2_quals, s2_reverse.DeletionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s3_quals, s3_reverse.DeletionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s4_quals, s4_reverse.DeletionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s5_quals, s5_reverse.DeletionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s6_quals, s6_reverse.DeletionQV(Orientation::NATIVE, false, false).Fastq());

    EXPECT_EQ(s1_quals, s1_reverse.InsertionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s2_quals, s2_reverse.InsertionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s3_quals, s3_reverse.InsertionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s4_quals, s4_reverse.InsertionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s5_quals, s5_reverse.InsertionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s6_quals, s6_reverse.InsertionQV(Orientation::NATIVE, false, false).Fastq());

    EXPECT_EQ(s1_quals, s1_reverse.MergeQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s2_quals, s2_reverse.MergeQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s3_quals, s3_reverse.MergeQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s4_quals, s4_reverse.MergeQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s5_quals, s5_reverse.MergeQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s6_quals, s6_reverse.MergeQV(Orientation::NATIVE, false, false).Fastq());

    EXPECT_EQ(s1_quals, s1_reverse.SubstitutionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s2_quals, s2_reverse.SubstitutionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s3_quals, s3_reverse.SubstitutionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s4_quals, s4_reverse.SubstitutionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s5_quals, s5_reverse.SubstitutionQV(Orientation::NATIVE, false, false).Fastq());
    EXPECT_EQ(s6_quals, s6_reverse.SubstitutionQV(Orientation::NATIVE, false, false).Fastq());

    //  - aligned, unclipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revquals_aligned, s1_reverse.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_revquals_aligned, s2_reverse.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_revquals_aligned, s3_reverse.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_revquals_aligned, s4_reverse.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_revquals_aligned, s5_reverse.DeletionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_revquals_aligned, s6_reverse.DeletionQV(Orientation::GENOMIC, true, false).Fastq());

    EXPECT_EQ(s1_revquals_aligned, s1_reverse.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_revquals_aligned, s2_reverse.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_revquals_aligned, s3_reverse.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_revquals_aligned, s4_reverse.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_revquals_aligned, s5_reverse.InsertionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_revquals_aligned, s6_reverse.InsertionQV(Orientation::GENOMIC, true, false).Fastq());

    EXPECT_EQ(s1_revquals_aligned, s1_reverse.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_revquals_aligned, s2_reverse.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_revquals_aligned, s3_reverse.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_revquals_aligned, s4_reverse.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_revquals_aligned, s5_reverse.MergeQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_revquals_aligned, s6_reverse.MergeQV(Orientation::GENOMIC, true, false).Fastq());

    EXPECT_EQ(s1_revquals_aligned, s1_reverse.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s2_revquals_aligned, s2_reverse.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s3_revquals_aligned, s3_reverse.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s4_revquals_aligned, s4_reverse.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s5_revquals_aligned, s5_reverse.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());
    EXPECT_EQ(s6_revquals_aligned, s6_reverse.SubstitutionQV(Orientation::GENOMIC, true, false).Fastq());

    //  - aligned, unclipped native output (native input)
    EXPECT_EQ(s1_quals_aligned, s1_reverse.DeletionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2_reverse.DeletionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3_reverse.DeletionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4_reverse.DeletionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5_reverse.DeletionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6_reverse.DeletionQV(Orientation::NATIVE, true, false).Fastq());

    EXPECT_EQ(s1_quals_aligned, s1_reverse.InsertionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2_reverse.InsertionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3_reverse.InsertionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4_reverse.InsertionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5_reverse.InsertionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6_reverse.InsertionQV(Orientation::NATIVE, true, false).Fastq());

    EXPECT_EQ(s1_quals_aligned, s1_reverse.MergeQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2_reverse.MergeQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3_reverse.MergeQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4_reverse.MergeQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5_reverse.MergeQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6_reverse.MergeQV(Orientation::NATIVE, true, false).Fastq());

    EXPECT_EQ(s1_quals_aligned, s1_reverse.SubstitutionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s2_quals_aligned, s2_reverse.SubstitutionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s3_quals_aligned, s3_reverse.SubstitutionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s4_quals_aligned, s4_reverse.SubstitutionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s5_quals_aligned, s5_reverse.SubstitutionQV(Orientation::NATIVE, true, false).Fastq());
    EXPECT_EQ(s6_quals_aligned, s6_reverse.SubstitutionQV(Orientation::NATIVE, true, false).Fastq());

    //  - aligned, clipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revquals_aligned_clipped, s1_reverse.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_revquals_aligned_clipped, s2_reverse.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_revquals_aligned_clipped, s3_reverse.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_revquals_aligned_clipped, s4_reverse.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_revquals_aligned_clipped, s5_reverse.DeletionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_revquals_aligned_clipped, s6_reverse.DeletionQV(Orientation::GENOMIC, true, true).Fastq());

    EXPECT_EQ(s1_revquals_aligned_clipped, s1_reverse.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_revquals_aligned_clipped, s2_reverse.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_revquals_aligned_clipped, s3_reverse.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_revquals_aligned_clipped, s4_reverse.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_revquals_aligned_clipped, s5_reverse.InsertionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_revquals_aligned_clipped, s6_reverse.InsertionQV(Orientation::GENOMIC, true, true).Fastq());

    EXPECT_EQ(s1_revquals_aligned_clipped, s1_reverse.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_revquals_aligned_clipped, s2_reverse.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_revquals_aligned_clipped, s3_reverse.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_revquals_aligned_clipped, s4_reverse.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_revquals_aligned_clipped, s5_reverse.MergeQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_revquals_aligned_clipped, s6_reverse.MergeQV(Orientation::GENOMIC, true, true).Fastq());

    EXPECT_EQ(s1_revquals_aligned_clipped, s1_reverse.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s2_revquals_aligned_clipped, s2_reverse.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s3_revquals_aligned_clipped, s3_reverse.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s4_revquals_aligned_clipped, s4_reverse.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s5_revquals_aligned_clipped, s5_reverse.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());
    EXPECT_EQ(s6_revquals_aligned_clipped, s6_reverse.SubstitutionQV(Orientation::GENOMIC, true, true).Fastq());

    //  - aligned, clipped native output (native input)
    EXPECT_EQ(s1_quals_aligned_clipped, s1_reverse.DeletionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2_reverse.DeletionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3_reverse.DeletionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4_reverse.DeletionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5_reverse.DeletionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6_reverse.DeletionQV(Orientation::NATIVE, true, true).Fastq());

    EXPECT_EQ(s1_quals_aligned_clipped, s1_reverse.InsertionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2_reverse.InsertionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3_reverse.InsertionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4_reverse.InsertionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5_reverse.InsertionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6_reverse.InsertionQV(Orientation::NATIVE, true, true).Fastq());

    EXPECT_EQ(s1_quals_aligned_clipped, s1_reverse.MergeQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2_reverse.MergeQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3_reverse.MergeQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4_reverse.MergeQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5_reverse.MergeQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6_reverse.MergeQV(Orientation::NATIVE, true, true).Fastq());

    EXPECT_EQ(s1_quals_aligned_clipped, s1_reverse.SubstitutionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s2_quals_aligned_clipped, s2_reverse.SubstitutionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s3_quals_aligned_clipped, s3_reverse.SubstitutionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s4_quals_aligned_clipped, s4_reverse.SubstitutionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s5_quals_aligned_clipped, s5_reverse.SubstitutionQV(Orientation::NATIVE, true, true).Fastq());
    EXPECT_EQ(s6_quals_aligned_clipped, s6_reverse.SubstitutionQV(Orientation::NATIVE, true, true).Fastq());
}

TEST(BamRecordTest, BaseTagsClippedAndAligned)
{
    const string s1_cigar = "4M3D4M";
    const string s2_cigar = "4M1D2I2D4M";
    const string s3_cigar = "4M1D2P2I2P2D4M";
    const string s4_cigar = "3S4M3D4M3S";
    const string s5_cigar = "2H4M3D4M3H";
    const string s6_cigar = "2H3S4M3D4M3S3H";

    const string s1_seq = "AACCGTTA";
    const string s2_seq = "ATCCTAGGTT";
    const string s3_seq = "ATCCTAGGTT";
    const string s4_seq = "TTTAACCGTTACCG";
    const string s5_seq = "AACCGTTA";
    const string s6_seq = "TTTAACCGTTACCG";

    const string s1_seq_aligned = "AACC---GTTA";
    const string s2_seq_aligned = "ATCC-TA--GGTT";
    const string s3_seq_aligned = "ATCC-**TA**--GGTT";
    const string s4_seq_aligned = "TTTAACC---GTTACCG";
    const string s5_seq_aligned = "AACC---GTTA";
    const string s6_seq_aligned = "TTTAACC---GTTACCG";

    const string s1_seq_aligned_clipped = "AACC---GTTA";
    const string s2_seq_aligned_clipped = "ATCC-TA--GGTT";
    const string s3_seq_aligned_clipped = "ATCC-**TA**--GGTT";
    const string s4_seq_aligned_clipped = "AACC---GTTA";
    const string s5_seq_aligned_clipped = "AACC---GTTA";
    const string s6_seq_aligned_clipped = "AACC---GTTA";

    const string s1_revseq = "TAACGGTT";
    const string s2_revseq = "AACCTAGGAT";
    const string s3_revseq = "AACCTAGGAT";
    const string s4_revseq = "CGGTAACGGTTAAA";
    const string s5_revseq = "TAACGGTT";
    const string s6_revseq = "CGGTAACGGTTAAA";

    const string s1_revseq_aligned = "TAAC---GGTT";
    const string s2_revseq_aligned = "AACC-TA--GGAT";
    const string s3_revseq_aligned = "AACC-**TA**--GGAT";
    const string s4_revseq_aligned = "CGGTAAC---GGTTAAA";
    const string s5_revseq_aligned = "TAAC---GGTT";
    const string s6_revseq_aligned = "CGGTAAC---GGTTAAA";

    const string s1_revseq_aligned_clipped = "TAAC---GGTT";
    const string s2_revseq_aligned_clipped = "AACC-TA--GGAT";
    const string s3_revseq_aligned_clipped = "AACC-**TA**--GGAT";
    const string s4_revseq_aligned_clipped = "TAAC---GGTT";
    const string s5_revseq_aligned_clipped = "TAAC---GGTT";
    const string s6_revseq_aligned_clipped = "TAAC---GGTT";

    const BamRecord s1 = tests::MakeCigaredBaseRecord(s1_seq, s1_cigar, false);
    const BamRecord s2 = tests::MakeCigaredBaseRecord(s2_seq, s2_cigar, false);
    const BamRecord s3 = tests::MakeCigaredBaseRecord(s3_seq, s3_cigar, false);
    const BamRecord s4 = tests::MakeCigaredBaseRecord(s4_seq, s4_cigar, false);
    const BamRecord s5 = tests::MakeCigaredBaseRecord(s5_seq, s5_cigar, false);
    const BamRecord s6 = tests::MakeCigaredBaseRecord(s6_seq, s6_cigar, false);
    const BamRecord s1_reverse = tests::MakeCigaredBaseRecord(s1_seq, s1_cigar, true);
    const BamRecord s2_reverse = tests::MakeCigaredBaseRecord(s2_seq, s2_cigar, true);
    const BamRecord s3_reverse = tests::MakeCigaredBaseRecord(s3_seq, s3_cigar, true);
    const BamRecord s4_reverse = tests::MakeCigaredBaseRecord(s4_seq, s4_cigar, true);
    const BamRecord s5_reverse = tests::MakeCigaredBaseRecord(s5_seq, s5_cigar, true);
    const BamRecord s6_reverse = tests::MakeCigaredBaseRecord(s6_seq, s6_cigar, true);

    // ----------------
    // forward strand
    // ----------------

    //  - "native" == "genomic"
    EXPECT_EQ(s1.DeletionTag(Orientation::NATIVE), s1.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2.DeletionTag(Orientation::NATIVE), s2.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3.DeletionTag(Orientation::NATIVE), s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s4.DeletionTag(Orientation::NATIVE), s4.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s5.DeletionTag(Orientation::NATIVE), s5.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s6.DeletionTag(Orientation::NATIVE), s6.DeletionTag(Orientation::GENOMIC));

    EXPECT_EQ(s1.SubstitutionTag(Orientation::NATIVE), s1.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2.SubstitutionTag(Orientation::NATIVE), s2.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3.SubstitutionTag(Orientation::NATIVE), s3.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s4.SubstitutionTag(Orientation::NATIVE), s4.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s5.SubstitutionTag(Orientation::NATIVE), s5.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s6.SubstitutionTag(Orientation::NATIVE), s6.SubstitutionTag(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output == genomic input
    EXPECT_EQ(s1_seq, s1.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_seq, s2.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_seq, s3.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s4_seq, s4.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s5_seq, s5.DeletionTag(Orientation::GENOMIC));
    EXPECT_EQ(s6_seq, s6.DeletionTag(Orientation::GENOMIC));

    EXPECT_EQ(s1_seq, s1.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s2_seq, s2.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s3_seq, s3.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s4_seq, s4.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s5_seq, s5.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_EQ(s6_seq, s6.SubstitutionTag(Orientation::GENOMIC));

    //  - aligned, unclipped genomic output == aligned, unclipped genomic input
    EXPECT_EQ(s1_seq_aligned, s1.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s2_seq_aligned, s2.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s3_seq_aligned, s3.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s4_seq_aligned, s4.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s5_seq_aligned, s5.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s6_seq_aligned, s6.DeletionTag(Orientation::GENOMIC, true, false));

    EXPECT_EQ(s1_seq_aligned, s1.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s2_seq_aligned, s2.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s3_seq_aligned, s3.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s4_seq_aligned, s4.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s5_seq_aligned, s5.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s6_seq_aligned, s6.SubstitutionTag(Orientation::GENOMIC, true, false));

    //  - aligned, clipped genomic output == aligned, clipped genomic output
    EXPECT_EQ(s1_seq_aligned_clipped, s1.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s2_seq_aligned_clipped, s2.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s3_seq_aligned_clipped, s3.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s4_seq_aligned_clipped, s4.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s5_seq_aligned_clipped, s5.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s6_seq_aligned_clipped, s6.DeletionTag(Orientation::GENOMIC, true, true));

    EXPECT_EQ(s1_seq_aligned_clipped, s1.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s2_seq_aligned_clipped, s2.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s3_seq_aligned_clipped, s3.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s4_seq_aligned_clipped, s4.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s5_seq_aligned_clipped, s5.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s6_seq_aligned_clipped, s6.SubstitutionTag(Orientation::GENOMIC, true, true));

    // ----------------
    // reverse strand
    // ----------------

    //  - "native" != "genomic"
    EXPECT_NE(s1_reverse.DeletionTag(Orientation::NATIVE), s1_reverse.DeletionTag(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.DeletionTag(Orientation::NATIVE), s2_reverse.DeletionTag(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.DeletionTag(Orientation::NATIVE), s3_reverse.DeletionTag(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.DeletionTag(Orientation::NATIVE), s4_reverse.DeletionTag(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.DeletionTag(Orientation::NATIVE), s5_reverse.DeletionTag(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.DeletionTag(Orientation::NATIVE), s6_reverse.DeletionTag(Orientation::GENOMIC));

    EXPECT_NE(s1_reverse.SubstitutionTag(Orientation::NATIVE), s1_reverse.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.SubstitutionTag(Orientation::NATIVE), s2_reverse.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.SubstitutionTag(Orientation::NATIVE), s3_reverse.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.SubstitutionTag(Orientation::NATIVE), s4_reverse.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.SubstitutionTag(Orientation::NATIVE), s5_reverse.SubstitutionTag(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.SubstitutionTag(Orientation::NATIVE), s6_reverse.SubstitutionTag(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revseq, s1_reverse.DeletionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s2_revseq, s2_reverse.DeletionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s3_revseq, s3_reverse.DeletionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s4_revseq, s4_reverse.DeletionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s5_revseq, s5_reverse.DeletionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s6_revseq, s6_reverse.DeletionTag(Orientation::GENOMIC, false, false));

    EXPECT_EQ(s1_revseq, s1_reverse.SubstitutionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s2_revseq, s2_reverse.SubstitutionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s3_revseq, s3_reverse.SubstitutionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s4_revseq, s4_reverse.SubstitutionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s5_revseq, s5_reverse.SubstitutionTag(Orientation::GENOMIC, false, false));
    EXPECT_EQ(s6_revseq, s6_reverse.SubstitutionTag(Orientation::GENOMIC, false, false));

    //  - unaligned, unclipped native output (native input)
    EXPECT_EQ(s1_seq, s1_reverse.DeletionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s2_seq, s2_reverse.DeletionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s3_seq, s3_reverse.DeletionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s4_seq, s4_reverse.DeletionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s5_seq, s5_reverse.DeletionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s6_seq, s6_reverse.DeletionTag(Orientation::NATIVE, false, false));

    EXPECT_EQ(s1_seq, s1_reverse.SubstitutionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s2_seq, s2_reverse.SubstitutionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s3_seq, s3_reverse.SubstitutionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s4_seq, s4_reverse.SubstitutionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s5_seq, s5_reverse.SubstitutionTag(Orientation::NATIVE, false, false));
    EXPECT_EQ(s6_seq, s6_reverse.SubstitutionTag(Orientation::NATIVE, false, false));

    //  - aligned, unclipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revseq_aligned, s1_reverse.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s2_revseq_aligned, s2_reverse.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s3_revseq_aligned, s3_reverse.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s4_revseq_aligned, s4_reverse.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s5_revseq_aligned, s5_reverse.DeletionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s6_revseq_aligned, s6_reverse.DeletionTag(Orientation::GENOMIC, true, false));

    EXPECT_EQ(s1_revseq_aligned, s1_reverse.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s2_revseq_aligned, s2_reverse.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s3_revseq_aligned, s3_reverse.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s4_revseq_aligned, s4_reverse.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s5_revseq_aligned, s5_reverse.SubstitutionTag(Orientation::GENOMIC, true, false));
    EXPECT_EQ(s6_revseq_aligned, s6_reverse.SubstitutionTag(Orientation::GENOMIC, true, false));

    //  - aligned, unclipped native output (native input)
    EXPECT_EQ(s1_seq_aligned, s1_reverse.DeletionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s2_seq_aligned, s2_reverse.DeletionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s3_seq_aligned, s3_reverse.DeletionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s4_seq_aligned, s4_reverse.DeletionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s5_seq_aligned, s5_reverse.DeletionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s6_seq_aligned, s6_reverse.DeletionTag(Orientation::NATIVE, true, false));

    EXPECT_EQ(s1_seq_aligned, s1_reverse.SubstitutionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s2_seq_aligned, s2_reverse.SubstitutionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s3_seq_aligned, s3_reverse.SubstitutionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s4_seq_aligned, s4_reverse.SubstitutionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s5_seq_aligned, s5_reverse.SubstitutionTag(Orientation::NATIVE, true, false));
    EXPECT_EQ(s6_seq_aligned, s6_reverse.SubstitutionTag(Orientation::NATIVE, true, false));

    //  - aligned, clipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revseq_aligned_clipped, s1_reverse.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s2_revseq_aligned_clipped, s2_reverse.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s3_revseq_aligned_clipped, s3_reverse.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s4_revseq_aligned_clipped, s4_reverse.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s5_revseq_aligned_clipped, s5_reverse.DeletionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s6_revseq_aligned_clipped, s6_reverse.DeletionTag(Orientation::GENOMIC, true, true));

    EXPECT_EQ(s1_revseq_aligned_clipped, s1_reverse.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s2_revseq_aligned_clipped, s2_reverse.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s3_revseq_aligned_clipped, s3_reverse.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s4_revseq_aligned_clipped, s4_reverse.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s5_revseq_aligned_clipped, s5_reverse.SubstitutionTag(Orientation::GENOMIC, true, true));
    EXPECT_EQ(s6_revseq_aligned_clipped, s6_reverse.SubstitutionTag(Orientation::GENOMIC, true, true));

    //  - aligned, clipped native output (native input)
    EXPECT_EQ(s1_seq_aligned_clipped, s1_reverse.DeletionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s2_seq_aligned_clipped, s2_reverse.DeletionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s3_seq_aligned_clipped, s3_reverse.DeletionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s4_seq_aligned_clipped, s4_reverse.DeletionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s5_seq_aligned_clipped, s5_reverse.DeletionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s6_seq_aligned_clipped, s6_reverse.DeletionTag(Orientation::NATIVE, true, true));

    EXPECT_EQ(s1_seq_aligned_clipped, s1_reverse.SubstitutionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s2_seq_aligned_clipped, s2_reverse.SubstitutionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s3_seq_aligned_clipped, s3_reverse.SubstitutionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s4_seq_aligned_clipped, s4_reverse.SubstitutionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s5_seq_aligned_clipped, s5_reverse.SubstitutionTag(Orientation::NATIVE, true, true));
    EXPECT_EQ(s6_seq_aligned_clipped, s6_reverse.SubstitutionTag(Orientation::NATIVE, true, true));

}

TEST(BamRecordTest, FrameTagsClippedAndAligned)
{
    const string s1_cigar = "4M3D4M";
    const string s2_cigar = "4M1D2I2D4M";
    const string s3_cigar = "4M1D2P2I2P2D4M";
    const string s4_cigar = "3S4M3D4M3S";
    const string s5_cigar = "2H4M3D4M3H";
    const string s6_cigar = "2H3S4M3D4M3S3H";

    typedef vector<uint16_t> f_data;

    const f_data s1_frames = { 10, 20, 10, 20, 10, 20, 10, 30 };
    const f_data s2_frames = { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 };
    const f_data s3_frames = { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 };
    const f_data s4_frames = { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 };
    const f_data s5_frames = { 10, 20, 10, 20, 10, 20, 10, 30 };
    const f_data s6_frames = { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 };

    const f_data s1_frames_aligned = { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s2_frames_aligned = { 10, 20, 10, 20, 0, 80, 70, 0, 0, 10, 20, 10, 30 };
    const f_data s3_frames_aligned = { 10, 20, 10, 20, 0, 0, 0, 80, 70, 0, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s4_frames_aligned = { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 };
    const f_data s5_frames_aligned = { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s6_frames_aligned = { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 };

    const f_data s1_frames_aligned_clipped = { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s2_frames_aligned_clipped = { 10, 20, 10, 20, 0, 80, 70, 0, 0, 10, 20, 10, 30 };
    const f_data s3_frames_aligned_clipped = { 10, 20, 10, 20, 0, 0, 0, 80, 70, 0, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s4_frames_aligned_clipped = { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s5_frames_aligned_clipped = { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 };
    const f_data s6_frames_aligned_clipped = { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 };

    const f_data s1_revframes = { 30, 10, 20, 10, 20, 10, 20, 10 };
    const f_data s2_revframes = { 30, 10, 20, 10, 70, 80, 20, 10, 20, 10 };
    const f_data s3_revframes = { 30, 10, 20, 10, 70, 80, 20, 10, 20, 10 };
    const f_data s4_revframes = { 50, 50, 50, 30, 10, 20, 10, 20, 10, 20, 10, 40, 40, 40 };
    const f_data s5_revframes = { 30, 10, 20, 10, 20, 10, 20, 10 };
    const f_data s6_revframes = { 50, 50, 50, 30, 10, 20, 10, 20, 10, 20, 10, 40, 40, 40 };

    const f_data s1_revframes_aligned = { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s2_revframes_aligned = { 30, 10, 20, 10, 0, 70, 80, 0, 0, 20, 10, 20, 10 };
    const f_data s3_revframes_aligned = { 30, 10, 20, 10, 0, 0, 0, 70, 80, 0, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s4_revframes_aligned = { 50, 50, 50, 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10, 40, 40, 40 };
    const f_data s5_revframes_aligned = { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s6_revframes_aligned = { 50, 50, 50, 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10, 40, 40, 40 };

    const f_data s1_revframes_aligned_clipped = { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s2_revframes_aligned_clipped = { 30, 10, 20, 10, 0, 70, 80, 0, 0, 20, 10, 20, 10 };
    const f_data s3_revframes_aligned_clipped = { 30, 10, 20, 10, 0, 0, 0, 70, 80, 0, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s4_revframes_aligned_clipped = { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s5_revframes_aligned_clipped = { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 };
    const f_data s6_revframes_aligned_clipped = { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 };

    const BamRecord s1 = tests::MakeCigaredFrameRecord(s1_frames, s1_cigar, false);
    const BamRecord s2 = tests::MakeCigaredFrameRecord(s2_frames, s2_cigar, false);
    const BamRecord s3 = tests::MakeCigaredFrameRecord(s3_frames, s3_cigar, false);
    const BamRecord s4 = tests::MakeCigaredFrameRecord(s4_frames, s4_cigar, false);
    const BamRecord s5 = tests::MakeCigaredFrameRecord(s5_frames, s5_cigar, false);
    const BamRecord s6 = tests::MakeCigaredFrameRecord(s6_frames, s6_cigar, false);
    const BamRecord s1_reverse = tests::MakeCigaredFrameRecord(s1_frames, s1_cigar, true);
    const BamRecord s2_reverse = tests::MakeCigaredFrameRecord(s2_frames, s2_cigar, true);
    const BamRecord s3_reverse = tests::MakeCigaredFrameRecord(s3_frames, s3_cigar, true);
    const BamRecord s4_reverse = tests::MakeCigaredFrameRecord(s4_frames, s4_cigar, true);
    const BamRecord s5_reverse = tests::MakeCigaredFrameRecord(s5_frames, s5_cigar, true);
    const BamRecord s6_reverse = tests::MakeCigaredFrameRecord(s6_frames, s6_cigar, true);

    // ----------------
    // forward strand
    // ----------------

    //  - "native" == "genomic"
    EXPECT_EQ(s1.IPD(Orientation::NATIVE), s1.IPD(Orientation::GENOMIC));
    EXPECT_EQ(s2.IPD(Orientation::NATIVE), s2.IPD(Orientation::GENOMIC));
    EXPECT_EQ(s3.IPD(Orientation::NATIVE), s3.IPD(Orientation::GENOMIC));
    EXPECT_EQ(s4.IPD(Orientation::NATIVE), s4.IPD(Orientation::GENOMIC));
    EXPECT_EQ(s5.IPD(Orientation::NATIVE), s5.IPD(Orientation::GENOMIC));
    EXPECT_EQ(s6.IPD(Orientation::NATIVE), s6.IPD(Orientation::GENOMIC));

    EXPECT_EQ(s1.PulseWidth(Orientation::NATIVE), s1.PulseWidth(Orientation::GENOMIC));
    EXPECT_EQ(s2.PulseWidth(Orientation::NATIVE), s2.PulseWidth(Orientation::GENOMIC));
    EXPECT_EQ(s3.PulseWidth(Orientation::NATIVE), s3.PulseWidth(Orientation::GENOMIC));
    EXPECT_EQ(s4.PulseWidth(Orientation::NATIVE), s4.PulseWidth(Orientation::GENOMIC));
    EXPECT_EQ(s5.PulseWidth(Orientation::NATIVE), s5.PulseWidth(Orientation::GENOMIC));
    EXPECT_EQ(s6.PulseWidth(Orientation::NATIVE), s6.PulseWidth(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output == genomic input
    EXPECT_EQ(s1_frames, s1.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(s2_frames, s2.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(s3_frames, s3.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(s4_frames, s4.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(s5_frames, s5.IPD(Orientation::GENOMIC).Data());
    EXPECT_EQ(s6_frames, s6.IPD(Orientation::GENOMIC).Data());

    EXPECT_EQ(s1_frames, s1.PulseWidth(Orientation::GENOMIC).Data());
    EXPECT_EQ(s2_frames, s2.PulseWidth(Orientation::GENOMIC).Data());
    EXPECT_EQ(s3_frames, s3.PulseWidth(Orientation::GENOMIC).Data());
    EXPECT_EQ(s4_frames, s4.PulseWidth(Orientation::GENOMIC).Data());
    EXPECT_EQ(s5_frames, s5.PulseWidth(Orientation::GENOMIC).Data());
    EXPECT_EQ(s6_frames, s6.PulseWidth(Orientation::GENOMIC).Data());

    //  - aligned, unclipped genomic output == aligned, unclipped genomic input
    EXPECT_EQ(s1_frames_aligned, s1.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s2_frames_aligned, s2.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s3_frames_aligned, s3.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s4_frames_aligned, s4.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s5_frames_aligned, s5.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s6_frames_aligned, s6.IPD(Orientation::GENOMIC, true, false).Data());

    EXPECT_EQ(s1_frames_aligned, s1.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s2_frames_aligned, s2.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s3_frames_aligned, s3.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s4_frames_aligned, s4.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s5_frames_aligned, s5.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s6_frames_aligned, s6.PulseWidth(Orientation::GENOMIC, true, false).Data());

    //  - aligned, clipped genomic output == aligned, clipped genomic output
    EXPECT_EQ(s1_frames_aligned_clipped, s1.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s2_frames_aligned_clipped, s2.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s3_frames_aligned_clipped, s3.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s4_frames_aligned_clipped, s4.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s5_frames_aligned_clipped, s5.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s6_frames_aligned_clipped, s6.IPD(Orientation::GENOMIC, true, true).Data());

    EXPECT_EQ(s1_frames_aligned_clipped, s1.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s2_frames_aligned_clipped, s2.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s3_frames_aligned_clipped, s3.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s4_frames_aligned_clipped, s4.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s5_frames_aligned_clipped, s5.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s6_frames_aligned_clipped, s6.PulseWidth(Orientation::GENOMIC, true, true).Data());

    // ----------------
    // reverse strand
    // ----------------

    //  - "native" != "genomic"
    EXPECT_NE(s1_reverse.IPD(Orientation::NATIVE), s1_reverse.IPD(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.IPD(Orientation::NATIVE), s2_reverse.IPD(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.IPD(Orientation::NATIVE), s3_reverse.IPD(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.IPD(Orientation::NATIVE), s4_reverse.IPD(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.IPD(Orientation::NATIVE), s5_reverse.IPD(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.IPD(Orientation::NATIVE), s6_reverse.IPD(Orientation::GENOMIC));

    EXPECT_NE(s1_reverse.PulseWidth(Orientation::NATIVE), s1_reverse.PulseWidth(Orientation::GENOMIC));
    EXPECT_NE(s2_reverse.PulseWidth(Orientation::NATIVE), s2_reverse.PulseWidth(Orientation::GENOMIC));
    EXPECT_NE(s3_reverse.PulseWidth(Orientation::NATIVE), s3_reverse.PulseWidth(Orientation::GENOMIC));
    EXPECT_NE(s4_reverse.PulseWidth(Orientation::NATIVE), s4_reverse.PulseWidth(Orientation::GENOMIC));
    EXPECT_NE(s5_reverse.PulseWidth(Orientation::NATIVE), s5_reverse.PulseWidth(Orientation::GENOMIC));
    EXPECT_NE(s6_reverse.PulseWidth(Orientation::NATIVE), s6_reverse.PulseWidth(Orientation::GENOMIC));

    //  - unaligned, unclipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revframes, s1_reverse.IPD(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s2_revframes, s2_reverse.IPD(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s3_revframes, s3_reverse.IPD(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s4_revframes, s4_reverse.IPD(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s5_revframes, s5_reverse.IPD(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s6_revframes, s6_reverse.IPD(Orientation::GENOMIC, false, false).Data());

    EXPECT_EQ(s1_revframes, s1_reverse.PulseWidth(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s2_revframes, s2_reverse.PulseWidth(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s3_revframes, s3_reverse.PulseWidth(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s4_revframes, s4_reverse.PulseWidth(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s5_revframes, s5_reverse.PulseWidth(Orientation::GENOMIC, false, false).Data());
    EXPECT_EQ(s6_revframes, s6_reverse.PulseWidth(Orientation::GENOMIC, false, false).Data());

    //  - unaligned, unclipped native output (native input)
    EXPECT_EQ(s1_frames, s1_reverse.IPD(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s2_frames, s2_reverse.IPD(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s3_frames, s3_reverse.IPD(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s4_frames, s4_reverse.IPD(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s5_frames, s5_reverse.IPD(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s6_frames, s6_reverse.IPD(Orientation::NATIVE, false, false).Data());

    EXPECT_EQ(s1_frames, s1_reverse.PulseWidth(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s2_frames, s2_reverse.PulseWidth(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s3_frames, s3_reverse.PulseWidth(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s4_frames, s4_reverse.PulseWidth(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s5_frames, s5_reverse.PulseWidth(Orientation::NATIVE, false, false).Data());
    EXPECT_EQ(s6_frames, s6_reverse.PulseWidth(Orientation::NATIVE, false, false).Data());

    //  - aligned, unclipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revframes_aligned, s1_reverse.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s2_revframes_aligned, s2_reverse.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s3_revframes_aligned, s3_reverse.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s4_revframes_aligned, s4_reverse.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s5_revframes_aligned, s5_reverse.IPD(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s6_revframes_aligned, s6_reverse.IPD(Orientation::GENOMIC, true, false).Data());

    EXPECT_EQ(s1_revframes_aligned, s1_reverse.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s2_revframes_aligned, s2_reverse.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s3_revframes_aligned, s3_reverse.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s4_revframes_aligned, s4_reverse.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s5_revframes_aligned, s5_reverse.PulseWidth(Orientation::GENOMIC, true, false).Data());
    EXPECT_EQ(s6_revframes_aligned, s6_reverse.PulseWidth(Orientation::GENOMIC, true, false).Data());

    //  - aligned, unclipped native output (native input)
    EXPECT_EQ(s1_frames_aligned, s1_reverse.IPD(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s2_frames_aligned, s2_reverse.IPD(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s3_frames_aligned, s3_reverse.IPD(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s4_frames_aligned, s4_reverse.IPD(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s5_frames_aligned, s5_reverse.IPD(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s6_frames_aligned, s6_reverse.IPD(Orientation::NATIVE, true, false).Data());

    EXPECT_EQ(s1_frames_aligned, s1_reverse.PulseWidth(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s2_frames_aligned, s2_reverse.PulseWidth(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s3_frames_aligned, s3_reverse.PulseWidth(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s4_frames_aligned, s4_reverse.PulseWidth(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s5_frames_aligned, s5_reverse.PulseWidth(Orientation::NATIVE, true, false).Data());
    EXPECT_EQ(s6_frames_aligned, s6_reverse.PulseWidth(Orientation::NATIVE, true, false).Data());

    //  - aligned, clipped genomic output (rev-comp of native input)
    EXPECT_EQ(s1_revframes_aligned_clipped, s1_reverse.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s2_revframes_aligned_clipped, s2_reverse.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s3_revframes_aligned_clipped, s3_reverse.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s4_revframes_aligned_clipped, s4_reverse.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s5_revframes_aligned_clipped, s5_reverse.IPD(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s6_revframes_aligned_clipped, s6_reverse.IPD(Orientation::GENOMIC, true, true).Data());

    EXPECT_EQ(s1_revframes_aligned_clipped, s1_reverse.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s2_revframes_aligned_clipped, s2_reverse.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s3_revframes_aligned_clipped, s3_reverse.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s4_revframes_aligned_clipped, s4_reverse.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s5_revframes_aligned_clipped, s5_reverse.PulseWidth(Orientation::GENOMIC, true, true).Data());
    EXPECT_EQ(s6_revframes_aligned_clipped, s6_reverse.PulseWidth(Orientation::GENOMIC, true, true).Data());

    //  - aligned, clipped native output (native input)
    EXPECT_EQ(s1_frames_aligned_clipped, s1_reverse.IPD(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s2_frames_aligned_clipped, s2_reverse.IPD(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s3_frames_aligned_clipped, s3_reverse.IPD(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s4_frames_aligned_clipped, s4_reverse.IPD(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s5_frames_aligned_clipped, s5_reverse.IPD(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s6_frames_aligned_clipped, s6_reverse.IPD(Orientation::NATIVE, true, true).Data());

    EXPECT_EQ(s1_frames_aligned_clipped, s1_reverse.PulseWidth(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s2_frames_aligned_clipped, s2_reverse.PulseWidth(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s3_frames_aligned_clipped, s3_reverse.PulseWidth(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s4_frames_aligned_clipped, s4_reverse.PulseWidth(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s5_frames_aligned_clipped, s5_reverse.PulseWidth(Orientation::NATIVE, true, true).Data());
    EXPECT_EQ(s6_frames_aligned_clipped, s6_reverse.PulseWidth(Orientation::NATIVE, true, true).Data());
}

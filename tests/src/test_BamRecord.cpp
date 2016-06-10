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
#include <array>
#include <initializer_list>
#include <string>
#include <vector>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

namespace tests {

static
BamRecordImpl CreateBamImpl(void)
{
    TagCollection tags;
    tags["HX"] = string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = vector<uint8_t>({34, 5, 125});
    tags["XY"] = static_cast<int32_t>(-42);

    BamRecordImpl bam;
    bam.Bin(42);
    bam.Flag(42);
    bam.InsertSize(42);
    bam.MapQuality(42);
    bam.MatePosition(42);
    bam.MateReferenceId(42);
    bam.Position(42);
    bam.ReferenceId(42);
    bam.Tags(tags);
    return bam;
}

static inline
BamRecord CreateBam(void)
{ return BamRecord{ CreateBamImpl() }; }

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

static inline
void CheckRawData(const BamRecord& bam)
{ CheckRawData(bam.impl_); }

static
BamRecordImpl MakeCigaredImpl(const string& seq,
                              const string& cigar,
                              const Strand strand)
{
    BamRecordImpl impl;
    impl.SetMapped(true).ReferenceId(0).Position(0).MapQuality(0);
    impl.CigarData(Cigar::FromStdString(cigar));
    impl.MateReferenceId(-1).MatePosition(-1).InsertSize(0);
    impl.SetSequenceAndQualities(seq, string(seq.size(), '*'));
    impl.SetReverseStrand(strand == Strand::REVERSE);
    return impl;
}

static inline
BamRecord MakeCigaredRecord(const string& seq,
                            const string& cigar,
                            const Strand strand)
{ return BamRecord{ MakeCigaredImpl(seq, cigar, strand) }; }

static
BamRecord MakeCigaredBaseRecord(const string& bases,
                                const string& cigar,
                                const Strand strand)
{
    TagCollection tags;
    tags["dt"] = bases;
    tags["st"] = bases;

    const string seq = string(bases.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredFrameRecord(const vector<uint16_t>& frames,
                                 const string& cigar,
                                 const Strand strand)
{
    TagCollection tags;
    tags["ip"] = frames;
    tags["pw"] = frames;

    const string seq = string(frames.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredQualRecord(const string& quals,
                                const string& cigar,
                                const Strand strand)
{
    TagCollection tags;
    tags["dq"] = quals;
    tags["iq"] = quals;
    tags["mq"] = quals;
    tags["sq"] = quals;

    const string seq = string(quals.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

// ----------------------------------------------------------
// helper structs and methods for checking combinations of:
//   aligned strand, orientation requested, alignment, clipping
// ----------------------------------------------------------

// generic result holder for various requested states
template<typename T>
struct ExpectedResult
{
public:
    ExpectedResult(std::initializer_list<T> init)
        : d_(init)
    {
        assert(12 == init.size());
    }

    T ForwardGenomic(void) const               { return d_.at(0); }
    T ForwardNative(void) const                { return d_.at(1); }
    T ForwardGenomicAligned(void) const        { return d_.at(2); }
    T ForwardNativeAligned(void) const         { return d_.at(3); }
    T ForwardGenomicAlignedClipped(void) const { return d_.at(4); }
    T ForwardNativeAlignedClipped(void) const  { return d_.at(5); }
    T ReverseGenomic(void) const               { return d_.at(6); }
    T ReverseNative(void) const                { return d_.at(7); }
    T ReverseGenomicAligned(void) const        { return d_.at(8); }
    T ReverseNativeAligned(void) const         { return d_.at(9); }
    T ReverseGenomicAlignedClipped(void) const { return d_.at(10); }
    T ReverseNativeAlignedClipped(void) const  { return d_.at(11); }

private:
    vector<T> d_;
};

// generic data type checker on the various requested states
template<typename DataType, typename MakeRecordType, typename FetchDataType>
void CheckAlignAndClip(const string& cigar,
                       const DataType& input,
                       const tests::ExpectedResult<DataType>& e,
                       const MakeRecordType& makeRecord,
                       const FetchDataType& fetchData)
{
    {   // map to forward strand
        const BamRecord b = makeRecord(input, cigar, Strand::FORWARD);
        EXPECT_EQ(e.ForwardGenomic(),               fetchData(b, Orientation::GENOMIC, false, false));
        EXPECT_EQ(e.ForwardNative(),                fetchData(b, Orientation::NATIVE,  false, false));
        EXPECT_EQ(e.ForwardGenomicAligned(),        fetchData(b, Orientation::GENOMIC, true,  false));
        EXPECT_EQ(e.ForwardNativeAligned(),         fetchData(b, Orientation::NATIVE,  true,  false));
        EXPECT_EQ(e.ForwardGenomicAlignedClipped(), fetchData(b, Orientation::GENOMIC, true,  true));
        EXPECT_EQ(e.ForwardNativeAlignedClipped(),  fetchData(b, Orientation::NATIVE,  true,  true));
    }
    {   // map to reverse strand
        const BamRecord b = makeRecord(input, cigar, Strand::REVERSE);
        EXPECT_EQ(e.ReverseGenomic(),               fetchData(b, Orientation::GENOMIC, false, false));
        EXPECT_EQ(e.ReverseNative(),                fetchData(b, Orientation::NATIVE,  false, false));
        EXPECT_EQ(e.ReverseGenomicAligned(),        fetchData(b, Orientation::GENOMIC, true,  false));
        EXPECT_EQ(e.ReverseNativeAligned(),         fetchData(b, Orientation::NATIVE,  true,  false));
        EXPECT_EQ(e.ReverseGenomicAlignedClipped(), fetchData(b, Orientation::GENOMIC, true,  true));
        EXPECT_EQ(e.ReverseNativeAlignedClipped(),  fetchData(b, Orientation::NATIVE,  true,  true));
    }
}

static
void CheckBaseTagsClippedAndAligned(const string& cigar,
                                    const string& input,
                                    const ExpectedResult<string>& e)
{
    // aligned record + DeletionTag, SubstitutionTag
    auto makeRecord = [](const string& bases,
                         const string& cigar,
                         const Strand strand)
    { return MakeCigaredBaseRecord(bases, cigar, strand); };

    // DeletionTag
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.DeletionTag(orientation, aligned, exciseSoftClips); }
    );

    // SubstitutionTag
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.SubstitutionTag(orientation, aligned, exciseSoftClips); }
    );
}

static
void CheckFrameTagsClippedAndAligned(const string& cigar,
                                     const vector<uint16_t>& input,
                                     const ExpectedResult<vector<uint16_t> >& e)
{

    // aligned record + IPD, PulseWidth
    auto makeRecord = [](const vector<uint16_t>& frames,
                         const string& cigar,
                         const Strand strand)
    { return tests::MakeCigaredFrameRecord(frames, cigar, strand); };

    // IPD
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.IPD(orientation, aligned, exciseSoftClips).Data(); }
    );

    // PulseWidth
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.PulseWidth(orientation, aligned, exciseSoftClips).Data(); }
    );
}

static
void CheckQualityTagsClippedAndAligned(const string& cigar,
                                       const string& input,
                                       const ExpectedResult<string>& e)
{
    // aligned record + DeletionQV, InsertionQV, MergeQV, SubstitutionQV
    auto makeRecord = [](const string& quals,
                         const string& cigar,
                         const Strand strand)
    { return tests::MakeCigaredQualRecord(quals, cigar, strand); };

    // DeletionQV
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.DeletionQV(orientation, aligned, exciseSoftClips).Fastq(); }
    );

    // InsertionQV
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.InsertionQV(orientation, aligned, exciseSoftClips).Fastq(); }
    );

    // MergeQV
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.MergeQV(orientation, aligned, exciseSoftClips).Fastq(); }
    );

    // SubstitutionQV
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.SubstitutionQV(orientation, aligned, exciseSoftClips).Fastq(); }
    );
}

static
void CheckQualitiesClippedAndAligned(const string& cigar,
                                     const string& input,
                                     const ExpectedResult<string>& e)
{
    // aligned record w/ dummy SEQ & QUALs under test
    auto makeRecord = [](const string& quals,
                         const string& cigar,
                         const Strand strand)
    {
        const string seq = string(quals.size(), 'N');
        auto record = tests::MakeCigaredRecord(seq, cigar, strand);
        record.Impl().SetSequenceAndQualities(seq, quals);
        return record;
    };

    // QUAL
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.Qualities(orientation, aligned, exciseSoftClips).Fastq(); }
    );
}

static
void CheckSequenceClippedAndAligned(const string& cigar,
                                    const string& input,
                                    const ExpectedResult<string>& e)
{
    // aligned record w/ SEQ
    auto makeRecord = [](const string& seq,
                         const string& cigar,
                         const Strand strand)
    { return tests::MakeCigaredRecord(seq, cigar, strand); };

    // SEQ
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.Sequence(orientation, aligned, exciseSoftClips); }
    );
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
    EXPECT_EQ(255, bam.impl_.MapQuality());
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
    EXPECT_THROW(bam.HoleNumber(), std::exception);
    EXPECT_FALSE(bam.HasNumPasses());
    EXPECT_THROW(bam.NumPasses(), std::exception);

    // 8888888888888888888888888888888888888
//    EXPECT_EQ(-1, bam.NumPasses());
//    EXPECT_EQ(-1, bam.QueryStart());
//    EXPECT_EQ(-1, bam.QueryEnd());
//    EXPECT_EQ(0, bam.ReadAccuracy());
//    EXPECT_EQ(-1, bam.ReferenceStart());
//    EXPECT_EQ(-1, bam.ReferenceEnd());
    // 8888888888888888888888888888888888888

    EXPECT_THROW(bam.DeletionQV(), std::exception);
    EXPECT_THROW(bam.DeletionTag(), std::exception);
    EXPECT_THROW(bam.InsertionQV(), std::exception);
    EXPECT_THROW(bam.MergeQV(), std::exception);
    EXPECT_THROW(bam.SubstitutionQV(), std::exception);
    EXPECT_THROW(bam.SubstitutionTag(), std::exception);

    // 8888888888888888888888888888888888888
//    EXPECT_FALSE(bam.HasDeletionQV());
//    EXPECT_FALSE(bam.HasDeletionTag());
//    EXPECT_FALSE(bam.HasInsertionQV());
//    EXPECT_FALSE(bam.HasMergeQV());
//    EXPECT_FALSE(bam.HasSubstitutionQV());
//    EXPECT_FALSE(bam.HasSubstitutionTag());

//    EXPECT_EQ(emptyString, bam.MovieName());
//    EXPECT_EQ(emptyString, bam.ReadGroupId());
    // 8888888888888888888888888888888888888

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
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    BamRecord bam3(move(tests::CreateBamImpl()));
#ifdef __clang__
#pragma clang diagnostic pop
#endif

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
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    bam4 = move(tests::CreateBamImpl());
#ifdef __clang__
#pragma clang diagnostic pop
#endif

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
    {
        SCOPED_TRACE("Simple CIGAR Sequence");
        tests::CheckSequenceClippedAndAligned(
            "13=",                  // CIGAR
            "ATATATCCCGGCG",        // input
            {
                "ATATATCCCGGCG",    // forward strand, genomic
                "ATATATCCCGGCG",    // forward strand, native
                "ATATATCCCGGCG",    // forward strand, genomic, aligned
                "ATATATCCCGGCG",    // forward strand, native,  aligned
                "ATATATCCCGGCG",    // forward strand, genomic, aligned + clipped
                "ATATATCCCGGCG",    // forward strand, native,  aligned + clipped
                "ATATATCCCGGCG",    // reverse strand, genomic
                "CGCCGGGATATAT",    // reverse strand, native
                "ATATATCCCGGCG",    // reverse strand, genomic, aligned
                "CGCCGGGATATAT",    // reverse strand, native,  aligned
                "ATATATCCCGGCG",    // reverse strand, genomic, aligned + clipped
                "CGCCGGGATATAT"     // reverse strand, native,  aligned + clipped
            }
        );
    }
}

TEST(BamRecordTest, QualitiesOrientation)
{
    {
        SCOPED_TRACE("Simple CIGAR Qualities");
        tests::CheckQualitiesClippedAndAligned(
            "13=",                  // CIGAR
            "?]?]?]?]?]?]*",        // input
            {
                "?]?]?]?]?]?]*",    // forward strand, genomic
                "?]?]?]?]?]?]*",    // forward strand, native
                "?]?]?]?]?]?]*",    // forward strand, genomic, aligned
                "?]?]?]?]?]?]*",    // forward strand, native,  aligned
                "?]?]?]?]?]?]*",    // forward strand, genomic, aligned + clipped
                "?]?]?]?]?]?]*",    // forward strand, native,  aligned + clipped
                "?]?]?]?]?]?]*",    // reverse strand, genomic
                "*]?]?]?]?]?]?",    // reverse strand, native
                "?]?]?]?]?]?]*",    // reverse strand, genomic, aligned
                "*]?]?]?]?]?]?",    // reverse strand, native,  aligned
                "?]?]?]?]?]?]*",    // reverse strand, genomic, aligned + clipped
                "*]?]?]?]?]?]?"     // reverse strand, native,  aligned + clipped
            }
        );
    }
}

TEST(BamRecordTest, SequenceTagsOrientation)
{
    {
        SCOPED_TRACE("Simple CIGAR Base Tags");
        tests::CheckBaseTagsClippedAndAligned(
            "13=",                  // CIGAR
            "ATATATCCCGGCG",        // input
            {
                "ATATATCCCGGCG",    // forward strand, genomic
                "ATATATCCCGGCG",    // forward strand, native
                "ATATATCCCGGCG",    // forward strand, genomic, aligned
                "ATATATCCCGGCG",    // forward strand, native, aligned
                "ATATATCCCGGCG",    // forward strand, genomic, aligned, clipped
                "ATATATCCCGGCG",    // forward strand, native, aligned, clipped
                "CGCCGGGATATAT",    // reverse strand, genomic
                "ATATATCCCGGCG",    // reverse strand, native
                "CGCCGGGATATAT",    // reverse strand, genomic, aligned
                "ATATATCCCGGCG",    // reverse strand, native, aligned
                "CGCCGGGATATAT",    // reverse strand, genomic, aligned, clipped
                "ATATATCCCGGCG"     // reverse strand, native, aligned, clipped
            }
        );
    }
}

TEST(BamRecordTest, FrameTagsOrientation)
{
    {
        SCOPED_TRACE("Simple CIGAR Frames");
        tests::CheckFrameTagsClippedAndAligned(
            "5=",                   // CIGAR
            { 0, 1, 2, 3, 4 },      // input
            {
                { 0, 1, 2, 3, 4 },  // forward strand, genomic
                { 0, 1, 2, 3, 4 },  // forward strand, native
                { 0, 1, 2, 3, 4 },  // forward strand, genomic, aligned
                { 0, 1, 2, 3, 4 },  // forward strand, native, aligned
                { 0, 1, 2, 3, 4 },  // forward strand, genomic, aligned, clipped
                { 0, 1, 2, 3, 4 },  // forward strand, native, aligned, clipped
                { 4, 3, 2, 1, 0 },  // reverse strand, genomic
                { 0, 1, 2, 3, 4 },  // reverse strand, native
                { 4, 3, 2, 1, 0 },  // reverse strand, genomic, aligned
                { 0, 1, 2, 3, 4 },  // reverse strand, native, aligned
                { 4, 3, 2, 1, 0 },  // reverse strand, genomic, aligned, clipped
                { 0, 1, 2, 3, 4 }   // reverse strand, native, aligned, clipped
            }
        );
    }
}

TEST(BamRecordTest, QualityTagsOrientation)
{
    {
        SCOPED_TRACE("Simple CIGAR Quality Tags");
        tests::CheckQualityTagsClippedAndAligned(
            "13=",                  // CIGAR
            "?]?]?]?]?]?]*",        // input
            {
                "?]?]?]?]?]?]*",    // forward strand, genomic
                "?]?]?]?]?]?]*",    // forward strand, native
                "?]?]?]?]?]?]*",    // forward strand, genomic, aligned
                "?]?]?]?]?]?]*",    // forward strand, native,  aligned
                "?]?]?]?]?]?]*",    // forward strand, genomic, aligned + clipped
                "?]?]?]?]?]?]*",    // forward strand, native,  aligned + clipped
                "*]?]?]?]?]?]?",    // reverse strand, genomic
                "?]?]?]?]?]?]*",    // reverse strand, native
                "*]?]?]?]?]?]?",    // reverse strand, genomic, aligned
                "?]?]?]?]?]?]*",    // reverse strand, native,  aligned
                "*]?]?]?]?]?]?",    // reverse strand, genomic, aligned + clipped
                "?]?]?]?]?]?]*"     // reverse strand, native,  aligned + clipped
            }
        );
    }
}

TEST(BamRecordTest, SequenceClippedAndAligned)
{
    {
        SCOPED_TRACE("CIGAR: 10=");
        tests::CheckSequenceClippedAndAligned(
            "10=",              // CIGAR
            "ATCCGCGGTT",       // input
            {
                "ATCCGCGGTT",   // forward strand, genomic
                "ATCCGCGGTT",   // forward strand, native
                "ATCCGCGGTT",   // forward strand, genomic, aligned
                "ATCCGCGGTT",   // forward strand, native,  aligned
                "ATCCGCGGTT",   // forward strand, genomic, aligned + clipped
                "ATCCGCGGTT",   // forward strand, native,  aligned + clipped
                "ATCCGCGGTT",   // reverse strand, genomic
                "AACCGCGGAT",   // reverse strand, native
                "ATCCGCGGTT",   // reverse strand, genomic, aligned
                "AACCGCGGAT",   // reverse strand, native,  aligned
                "ATCCGCGGTT",   // reverse strand, genomic, aligned + clipped
                "AACCGCGGAT"    // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 3=4N3=");
        tests::CheckSequenceClippedAndAligned(
            "3=4N3=",       // CIGAR
            "ACGTT",        // input
            {
                "ACGTT",    // forward strand, genomic
                "ACGTT",    // forward strand, native
                "ACGTT",    // forward strand, genomic, aligned
                "ACGTT",    // forward strand, native,  aligned
                "ACGTT",    // forward strand, genomic, aligned + clipped
                "ACGTT",    // forward strand, native,  aligned + clipped
                "ACGTT",    // reverse strand, genomic
                "AACGT",    // reverse strand, native
                "ACGTT",    // reverse strand, genomic, aligned
                "AACGT",    // reverse strand, native,  aligned
                "ACGTT",    // reverse strand, genomic, aligned + clipped
                "AACGT"     // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 1S8=1S");
        tests::CheckSequenceClippedAndAligned(
            "1S8=1S",           // CIGAR
            "ACCCGCGGTT",       // input
            {
                "ACCCGCGGTT",   // forward strand, genomic
                "ACCCGCGGTT",   // forward strand, native
                "ACCCGCGGTT",   // forward strand, genomic, aligned
                "ACCCGCGGTT",   // forward strand, native,  aligned
                "CCCGCGGT",     // forward strand, genomic, aligned + clipped
                "CCCGCGGT",     // forward strand, native,  aligned + clipped
                "ACCCGCGGTT",   // reverse strand, genomic
                "AACCGCGGGT",   // reverse strand, native
                "ACCCGCGGTT",   // reverse strand, genomic, aligned
                "AACCGCGGGT",   // reverse strand, native,  aligned
                "CCCGCGGT",     // reverse strand, genomic, aligned + clipped
                "ACCGCGGG"      // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 1H8=1H");
        tests::CheckSequenceClippedAndAligned(
            "1H8=1H",           // CIGAR
            "ATCGCGGT",         // input
            {
                "ATCGCGGT",     // forward strand, genomic
                "ATCGCGGT",     // forward strand, native
                "ATCGCGGT",     // forward strand, genomic, aligned
                "ATCGCGGT",     // forward strand, native,  aligned
                "ATCGCGGT",     // forward strand, genomic, aligned + clipped
                "ATCGCGGT",     // forward strand, native,  aligned + clipped
                "ATCGCGGT",     // reverse strand, genomic
                "ACCGCGAT",     // reverse strand, native
                "ATCGCGGT",     // reverse strand, genomic, aligned
                "ACCGCGAT",     // reverse strand, native,  aligned
                "ATCGCGGT",     // reverse strand, genomic, aligned + clipped
                "ACCGCGAT"      // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2S6=2S");
        tests::CheckSequenceClippedAndAligned(
            "2S6=2S",           // CIGAR
            "AGCCGCGGTT",       // input
            {
                "AGCCGCGGTT",   // forward strand, genomic
                "AGCCGCGGTT",   // forward strand, native
                "AGCCGCGGTT",   // forward strand, genomic, aligned
                "AGCCGCGGTT",   // forward strand, native,  aligned
                "CCGCGG",       // forward strand, genomic, aligned + clipped
                "CCGCGG",       // forward strand, native,  aligned + clipped
                "AGCCGCGGTT",   // reverse strand, genomic
                "AACCGCGGCT",   // reverse strand, native
                "AGCCGCGGTT",   // reverse strand, genomic, aligned
                "AACCGCGGCT",   // reverse strand, native,  aligned
                "CCGCGG",       // reverse strand, genomic, aligned + clipped
                "CCGCGG"        // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2S3=2I3=2S");
        tests::CheckSequenceClippedAndAligned(
            "2S3=2I3=2S",           // CIGAR
            "ATCCGNNCGGTT",         // input
            {
                "ATCCGNNCGGTT",     // forward strand, genomic
                "ATCCGNNCGGTT",     // forward strand, native
                "ATCCGNNCGGTT",     // forward strand, genomic, aligned
                "ATCCGNNCGGTT",     // forward strand, native,  aligned
                "CCGNNCGG",         // forward strand, genomic, aligned + clipped
                "CCGNNCGG",         // forward strand, native,  aligned + clipped
                "ATCCGNNCGGTT",     // reverse strand, genomic
                "AACCGNNCGGAT",     // reverse strand, native
                "ATCCGNNCGGTT",     // reverse strand, genomic, aligned
                "AACCGNNCGGAT",     // reverse strand, native,  aligned
                "CCGNNCGG",         // reverse strand, genomic, aligned + clipped
                "CCGNNCGG"          // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H6=2H");
        tests::CheckSequenceClippedAndAligned(
            "2H6=2H",       // CIGAR
            "CAGCGG",       // input
            {
                "CAGCGG",   // forward strand, genomic
                "CAGCGG",   // forward strand, native
                "CAGCGG",   // forward strand, genomic, aligned
                "CAGCGG",   // forward strand, native,  aligned
                "CAGCGG",   // forward strand, genomic, aligned + clipped
                "CAGCGG",   // forward strand, native,  aligned + clipped
                "CAGCGG",   // reverse strand, genomic
                "CCGCTG",   // reverse strand, native
                "CAGCGG",   // reverse strand, genomic, aligned
                "CCGCTG",   // reverse strand, native,  aligned
                "CAGCGG",   // reverse strand, genomic, aligned + clipped
                "CCGCTG"    // reverse strand, native,  aligned + clipped
            }
        );
    }
}

TEST(BamRecordTest, ClippingOrientationAndAlignment)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        tests::CheckSequenceClippedAndAligned(
            "4=3D4=",           // CIGAR
            "AACCGTTA",         // input
            {
                "AACCGTTA",     // forward strand, genomic
                "AACCGTTA",     // forward strand, native
                "AACC---GTTA",  // forward strand, genomic, aligned
                "AACC---GTTA",  // forward strand, native,  aligned
                "AACC---GTTA",  // forward strand, genomic, aligned + clipped
                "AACC---GTTA",  // forward strand, native,  aligned + clipped
                "AACCGTTA",     // reverse strand, genomic
                "TAACGGTT",     // reverse strand, native
                "AACC---GTTA",  // reverse strand, genomic, aligned
                "TAAC---GGTT",  // reverse strand, native,  aligned
                "AACC---GTTA",  // reverse strand, genomic, aligned + clipped
                "TAAC---GGTT"   // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2I2D4=");
        tests::CheckSequenceClippedAndAligned(
            "4=1D2I2D4=",           // CIGAR
            "ATCCTAGGTT",           // input
            {
                "ATCCTAGGTT",       // forward strand, genomic
                "ATCCTAGGTT",       // forward strand, native
                "ATCC-TA--GGTT",    // forward strand, genomic, aligned
                "ATCC-TA--GGTT",    // forward strand, native,  aligned
                "ATCC-TA--GGTT",    // forward strand, genomic, aligned + clipped
                "ATCC-TA--GGTT",    // forward strand, native,  aligned + clipped
                "ATCCTAGGTT",       // reverse strand, genomic
                "AACCTAGGAT",       // reverse strand, native
                "ATCC-TA--GGTT",    // reverse strand, genomic, aligned
                "AACC--TA-GGAT",    // reverse strand, native,  aligned
                "ATCC-TA--GGTT",    // reverse strand, genomic, aligned + clipped
                "AACC--TA-GGAT"     // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2P2I2P2D4=");
        tests::CheckSequenceClippedAndAligned(
            "4=1D2P2I2P2D4=",           // CIGAR
            "ATCCTAGGTT",               // input
            {
                "ATCCTAGGTT",           // forward strand, genomic
                "ATCCTAGGTT",           // forward strand, native
                "ATCC-**TA**--GGTT",    // forward strand, genomic, aligned
                "ATCC-**TA**--GGTT",    // forward strand, native,  aligned
                "ATCC-**TA**--GGTT",    // forward strand, genomic, aligned + clipped
                "ATCC-**TA**--GGTT",    // forward strand, native,  aligned + clipped
                "ATCCTAGGTT",           // reverse strand, genomic
                "AACCTAGGAT",           // reverse strand, native
                "ATCC-**TA**--GGTT",    // reverse strand, genomic, aligned
                "AACC--**TA**-GGAT",    // reverse strand, native,  aligned
                "ATCC-**TA**--GGTT",    // reverse strand, genomic, aligned + clipped
                "AACC--**TA**-GGAT"     // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2S4=3D4=3S");
        tests::CheckSequenceClippedAndAligned(
            "2S4=3D4=3S",               // CIGAR
            "TTAACCGTTACCG",            // input
            {
                "TTAACCGTTACCG",        // forward strand, genomic
                "TTAACCGTTACCG",        // forward strand, native
                "TTAACC---GTTACCG",     // forward strand, genomic, aligned
                "TTAACC---GTTACCG",     // forward strand, native,  aligned
                "AACC---GTTA",          // forward strand, genomic, aligned + clipped
                "AACC---GTTA",          // forward strand, native,  aligned + clipped
                "TTAACCGTTACCG",        // reverse strand, genomic
                "CGGTAACGGTTAA",        // reverse strand, native
                "TTAACC---GTTACCG",     // reverse strand, genomic, aligned
                "CGGTAAC---GGTTAA",     // reverse strand, native,  aligned
                "AACC---GTTA",          // reverse strand, genomic, aligned + clipped
                "TAAC---GGTT"           // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H4=3D4=3H");
        tests::CheckSequenceClippedAndAligned(
            "2H4=3D4=3H",       // CIGAR
            "AACCGTTA",         // input
            {
                "AACCGTTA",     // forward strand, genomic
                "AACCGTTA",     // forward strand, native
                "AACC---GTTA",  // forward strand, genomic, aligned
                "AACC---GTTA",  // forward strand, native,  aligned
                "AACC---GTTA",  // forward strand, genomic, aligned + clipped
                "AACC---GTTA",  // forward strand, native,  aligned + clipped
                "AACCGTTA",     // reverse strand, genomic
                "TAACGGTT",     // reverse strand, native
                "AACC---GTTA",  // reverse strand, genomic, aligned
                "TAAC---GGTT",  // reverse strand, native,  aligned
                "AACC---GTTA",  // reverse strand, genomic, aligned + clipped
                "TAAC---GGTT"   // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H2S4=3D4=3S3H");
        tests::CheckSequenceClippedAndAligned(
            "2H2S4=3D4=3S3H",           // CIGAR
            "TTAACCGTTACCG",            // input
            {
                "TTAACCGTTACCG",        // forward strand, genomic
                "TTAACCGTTACCG",        // forward strand, native
                "TTAACC---GTTACCG",     // forward strand, genomic, aligned
                "TTAACC---GTTACCG",     // forward strand, native,  aligned
                "AACC---GTTA",          // forward strand, genomic, aligned + clipped
                "AACC---GTTA",          // forward strand, native,  aligned + clipped
                "TTAACCGTTACCG",        // reverse strand, genomic
                "CGGTAACGGTTAA",        // reverse strand, native
                "TTAACC---GTTACCG",     // reverse strand, genomic, aligned
                "CGGTAAC---GGTTAA",     // reverse strand, native,  aligned
                "AACC---GTTA",          // reverse strand, genomic, aligned + clipped
                "TAAC---GGTT"           // reverse strand, native,  aligned + clipped
            }
        );
    }
}

TEST(BamRecordTest, QualityTagsClippedAndAligned)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        tests::CheckQualityTagsClippedAndAligned(
            "4=3D4=",           // CIGAR
            "?]?]?]?@",         // input
            {
                "?]?]?]?@",     // forward strand, genomic
                "?]?]?]?@",     // forward strand, native
                "?]?]!!!?]?@",  // forward strand, genomic, aligned
                "?]?]!!!?]?@",  // forward strand, native,  aligned
                "?]?]!!!?]?@",  // forward strand, genomic, aligned + clipped
                "?]?]!!!?]?@",  // forward strand, native,  aligned + clipped
                "@?]?]?]?",     // reverse strand, genomic
                "?]?]?]?@",     // reverse strand, native
                "@?]?!!!]?]?",  // reverse strand, genomic, aligned
                "?]?]!!!?]?@",  // reverse strand, native,  aligned
                "@?]?!!!]?]?",  // reverse strand, genomic, aligned + clipped
                "?]?]!!!?]?@"   // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2I2D4=");
        tests::CheckQualityTagsClippedAndAligned(
            "4=1D2I2D4=",           // CIGAR
            "?]?]87?]?@",           // input
            {
                "?]?]87?]?@",       // forward strand, genomic
                "?]?]87?]?@",       // forward strand, native
                "?]?]!87!!?]?@",    // forward strand, genomic, aligned
                "?]?]!87!!?]?@",    // forward strand, native,  aligned
                "?]?]!87!!?]?@",    // forward strand, genomic, aligned + clipped
                "?]?]!87!!?]?@",    // forward strand, native,  aligned + clipped
                "@?]?78]?]?",       // reverse strand, genomic
                "?]?]87?]?@",       // reverse strand, native
                "@?]?!78!!]?]?",    // reverse strand, genomic, aligned
                "?]?]!!87!?]?@",    // reverse strand, native,  aligned
                "@?]?!78!!]?]?",    // reverse strand, genomic, aligned + clipped
                "?]?]!!87!?]?@"     // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2P2I2P2D4=");
        tests::CheckQualityTagsClippedAndAligned(
            "4=1D2P2I2P2D4=",       // CIGAR
            "?]?]87?]?@",           // input
        {
            "?]?]87?]?@",           // forward strand, genomic
            "?]?]87?]?@",           // forward strand, native
            "?]?]!!!87!!!!?]?@",    // forward strand, genomic, aligned
            "?]?]!!!87!!!!?]?@",    // forward strand, native,  aligned
            "?]?]!!!87!!!!?]?@",    // forward strand, genomic, aligned + clipped
            "?]?]!!!87!!!!?]?@",    // forward strand, native,  aligned + clipped
            "@?]?78]?]?",           // reverse strand, genomic
            "?]?]87?]?@",           // reverse strand, native
            "@?]?!!!78!!!!]?]?",    // reverse strand, genomic, aligned
            "?]?]!!!!87!!!?]?@",    // reverse strand, native,  aligned
            "@?]?!!!78!!!!]?]?",    // reverse strand, genomic, aligned + clipped
            "?]?]!!!!87!!!?]?@"     // reverse strand, native,  aligned + clipped
        }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 3S4=3D4=3S");
        tests::CheckQualityTagsClippedAndAligned(
            "3S4=3D4=3S",               // CIGAR
            "vvv?]?]?]?@xxx",           // input
            {
                "vvv?]?]?]?@xxx",       // forward strand, genomic
                "vvv?]?]?]?@xxx",       // forward strand, native
                "vvv?]?]!!!?]?@xxx",    // forward strand, genomic, aligned
                "vvv?]?]!!!?]?@xxx",    // forward strand, native, aligned
                "?]?]!!!?]?@",          // forward strand, genomic, aligned, clipped
                "?]?]!!!?]?@",          // forward strand, native, aligned, clipped
                "xxx@?]?]?]?vvv",       // reverse strand, genomic
                "vvv?]?]?]?@xxx",       // reverse strand, native
                "xxx@?]?!!!]?]?vvv",    // reverse strand, genomic, aligned
                "vvv?]?]!!!?]?@xxx",    // reverse strand, native, aligned
                "@?]?!!!]?]?",          // reverse strand, genomic, aligned, clipped
                "?]?]!!!?]?@"           // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H4=3D4=3H");
        tests::CheckQualityTagsClippedAndAligned(
            "2H4=3D4=3H",       // CIGAR
            "?]?]?]?@",         // input
            {
                "?]?]?]?@",     // forward strand, genomic
                "?]?]?]?@",     // forward strand, native
                "?]?]!!!?]?@",  // forward strand, genomic, aligned
                "?]?]!!!?]?@",  // forward strand, native, aligned
                "?]?]!!!?]?@",  // forward strand, genomic, aligned, clipped
                "?]?]!!!?]?@",  // forward strand, native, aligned, clipped
                "@?]?]?]?",     // reverse strand, genomic
                "?]?]?]?@",     // reverse strand, native
                "@?]?!!!]?]?",  // reverse strand, genomic, aligned
                "?]?]!!!?]?@",  // reverse strand, native, aligned
                "@?]?!!!]?]?",  // reverse strand, genomic, aligned, clipped
                "?]?]!!!?]?@"   // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H3S4=3D4=3S3H");
        tests::CheckQualityTagsClippedAndAligned(
            "2H3S4=3D4=3S3H",           // CIGAR
            "vvv?]?]?]?@xxx",           // input
            {
                "vvv?]?]?]?@xxx",       // forward strand, genomic
                "vvv?]?]?]?@xxx",       // forward strand, native
                "vvv?]?]!!!?]?@xxx",    // forward strand, genomic, aligned
                "vvv?]?]!!!?]?@xxx",    // forward strand, native, aligned
                "?]?]!!!?]?@",          // forward strand, genomic, aligned, clipped
                "?]?]!!!?]?@",          // forward strand, native, aligned, clipped
                "xxx@?]?]?]?vvv",       // reverse strand, genomic
                "vvv?]?]?]?@xxx",       // reverse strand, native
                "xxx@?]?!!!]?]?vvv",    // reverse strand, genomic, aligned
                "vvv?]?]!!!?]?@xxx",    // reverse strand, native, aligned
                "@?]?!!!]?]?",          // reverse strand, genomic, aligned, clipped
                "?]?]!!!?]?@"           // reverse strand, native, aligned, clipped
            }
        );
    }
}

TEST(BamRecordTest, BaseTagsClippedAndAligned)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        tests::CheckBaseTagsClippedAndAligned(
            "4=3D4=",           // CIGAR
            "AACCGTTA",         // input
            {
                "AACCGTTA",     // forward strand, genomic
                "AACCGTTA",     // forward strand, native
                "AACC---GTTA",  // forward strand, genomic, aligned
                "AACC---GTTA",  // forward strand, native, aligned
                "AACC---GTTA",  // forward strand, genomic, aligned, clipped
                "AACC---GTTA",  // forward strand, native, aligned, clipped
                "TAACGGTT",     // reverse strand, genomic
                "AACCGTTA",     // reverse strand, native
                "TAAC---GGTT",  // reverse strand, genomic, aligned
                "AACC---GTTA",  // reverse strand, native, aligned
                "TAAC---GGTT",  // reverse strand, genomic, aligned, clipped
                "AACC---GTTA"   // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2I2D4=");
        tests::CheckBaseTagsClippedAndAligned(
            "4=1D2I2D4=",           // CIGAR
            "ATCCTAGGTT",           // input
            {
                "ATCCTAGGTT",       // forward strand, genomic
                "ATCCTAGGTT",       // forward strand, native
                "ATCC-TA--GGTT",    // forward strand, genomic, aligned
                "ATCC-TA--GGTT",    // forward strand, native, aligned
                "ATCC-TA--GGTT",    // forward strand, genomic, aligned, clipped
                "ATCC-TA--GGTT",    // forward strand, native, aligned, clipped
                "AACCTAGGAT",       // reverse strand, genomic
                "ATCCTAGGTT",       // reverse strand, native
                "AACC-TA--GGAT",    // reverse strand, genomic, aligned
                "ATCC--TA-GGTT",    // reverse strand, native, aligned
                "AACC-TA--GGAT",    // reverse strand, genomic, aligned, clipped
                "ATCC--TA-GGTT"     // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2P2I2P2D4=");
        tests::CheckBaseTagsClippedAndAligned(
            "4=1D2P2I2P2D4=",           // CIGAR
            "ATCCTAGGTT",               // input
            {
                "ATCCTAGGTT",           // forward strand, genomic
                "ATCCTAGGTT",           // forward strand, native
                "ATCC-**TA**--GGTT",    // forward strand, genomic, aligned
                "ATCC-**TA**--GGTT",    // forward strand, native, aligned
                "ATCC-**TA**--GGTT",    // forward strand, genomic, aligned, clipped
                "ATCC-**TA**--GGTT",    // forward strand, native, aligned, clipped
                "AACCTAGGAT",           // reverse strand, genomic
                "ATCCTAGGTT",           // reverse strand, native
                "AACC-**TA**--GGAT",    // reverse strand, genomic, aligned
                "ATCC--**TA**-GGTT",    // reverse strand, native, aligned
                "AACC-**TA**--GGAT",    // reverse strand, genomic, aligned, clipped
                "ATCC--**TA**-GGTT"     // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 3S4=3D4=3S");
        tests::CheckBaseTagsClippedAndAligned(
            "3S4=3D4=3S",               // CIGAR
            "TTTAACCGTTACCG",           // input
            {
                "TTTAACCGTTACCG",       // forward strand, genomic
                "TTTAACCGTTACCG",       // forward strand, native
                "TTTAACC---GTTACCG",    // forward strand, genomic, aligned
                "TTTAACC---GTTACCG",    // forward strand, native, aligned
                "AACC---GTTA",          // forward strand, genomic, aligned, clipped
                "AACC---GTTA",          // forward strand, native, aligned, clipped
                "CGGTAACGGTTAAA",       // reverse strand, genomic
                "TTTAACCGTTACCG",       // reverse strand, native
                "CGGTAAC---GGTTAAA",    // reverse strand, genomic, aligned
                "TTTAACC---GTTACCG",    // reverse strand, native, aligned
                "TAAC---GGTT",          // reverse strand, genomic, aligned, clipped
                "AACC---GTTA"           // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H4=3D4=3H");
        tests::CheckBaseTagsClippedAndAligned(
            "2H4=3D4=3H",       // CIGAR
            "AACCGTTA",         // input
            {
                "AACCGTTA",     // forward strand, genomic
                "AACCGTTA",     // forward strand, native
                "AACC---GTTA",  // forward strand, genomic, aligned
                "AACC---GTTA",  // forward strand, native, aligned
                "AACC---GTTA",  // forward strand, genomic, aligned, clipped
                "AACC---GTTA",  // forward strand, native, aligned, clipped
                "TAACGGTT",     // reverse strand, genomic
                "AACCGTTA",     // reverse strand, native
                "TAAC---GGTT",  // reverse strand, genomic, aligned
                "AACC---GTTA",  // reverse strand, native, aligned
                "TAAC---GGTT",  // reverse strand, genomic, aligned, clipped
                "AACC---GTTA"   // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H3S4=3D4=3S3H");
        tests::CheckBaseTagsClippedAndAligned(
            "2H3S4=3D4=3S3H",           // CIGAR
            "TTTAACCGTTACCG",           // input
            {
                "TTTAACCGTTACCG",       // forward strand, genomic
                "TTTAACCGTTACCG",       // forward strand, native
                "TTTAACC---GTTACCG",    // forward strand, genomic, aligned
                "TTTAACC---GTTACCG",    // forward strand, native, aligned
                "AACC---GTTA",          // forward strand, genomic, aligned, clipped
                "AACC---GTTA",          // forward strand, native, aligned, clipped
                "CGGTAACGGTTAAA",       // reverse strand, genomic
                "TTTAACCGTTACCG",       // reverse strand, native
                "CGGTAAC---GGTTAAA",    // reverse strand, genomic, aligned
                "TTTAACC---GTTACCG",    // reverse strand, native, aligned
                "TAAC---GGTT",          // reverse strand, genomic, aligned, clipped
                "AACC---GTTA"           // reverse strand, native, aligned, clipped
            }
        );
    }
}

TEST(BamRecordTest, FrameTagsClippedAndAligned)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        tests::CheckFrameTagsClippedAndAligned(
            "4=3D4=",                                           // CIGAR
            { 10, 20, 10, 20, 10, 20, 10, 30 },                 // input
            {
                { 10, 20, 10, 20, 10, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 10, 20, 10, 20, 10, 30 },             // forward strand, native
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 10, 20, 10, 20, 10 },             // reverse strand, genomic
                { 10, 20, 10, 20, 10, 20, 10, 30 },             // reverse strand, native
                { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // reverse strand, native, aligned
                { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 }     // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2I2D4=");
        tests::CheckFrameTagsClippedAndAligned(
            "4=1D2I2D4=",                                               // CIGAR
            { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },                 // input
            {
                { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },             // forward strand, native
                { 10, 20, 10, 20, 0, 80, 70, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned
                { 10, 20, 10, 20, 0, 80, 70, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned
                { 10, 20, 10, 20, 0, 80, 70, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 80, 70, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 10, 70, 80, 20, 10, 20, 10 },             // reverse strand, genomic
                { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },             // reverse strand, native
                { 30, 10, 20, 10, 0, 70, 80, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned
                { 10, 20, 10, 20, 0, 0, 80, 70, 0, 10, 20, 10, 30 },    // reverse strand, native, aligned
                { 30, 10, 20, 10, 0, 70, 80, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 80, 70, 0, 10, 20, 10, 30 }     // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 4=1D2P2I2P2D4=");
        tests::CheckFrameTagsClippedAndAligned(
            "4=1D2P2I2P2D4=",                                                   // CIGAR
            { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },                         // input
        {
            { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },                         // forward strand, genomic
            { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },                         // forward strand, native
            { 10, 20, 10, 20, 0, 0, 0, 80, 70, 0, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned
            { 10, 20, 10, 20, 0, 0, 0, 80, 70, 0, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned
            { 10, 20, 10, 20, 0, 0, 0, 80, 70, 0, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned, clipped
            { 10, 20, 10, 20, 0, 0, 0, 80, 70, 0, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned, clipped
            { 30, 10, 20, 10, 70, 80, 20, 10, 20, 10 },                         // reverse strand, genomic
            { 10, 20, 10, 20, 80, 70, 10, 20, 10, 30 },                         // reverse strand, native
            { 30, 10, 20, 10, 0, 0, 0, 70, 80, 0, 0, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned
            { 10, 20, 10, 20, 0, 0, 0, 0, 80, 70, 0, 0, 0, 10, 20, 10, 30 },    // reverse strand, native, aligned
            { 30, 10, 20, 10, 0, 0, 0, 70, 80, 0, 0, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned, clipped
            { 10, 20, 10, 20, 0, 0, 0, 0, 80, 70, 0, 0, 0, 10, 20, 10, 30 }     // reverse strand, native, aligned, clipped
        }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 3S4=3D4=3S");
        tests::CheckFrameTagsClippedAndAligned(
            "3S4=3D4=3S",                                                               // CIGAR
            { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },                 // input
            {
                { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },             // forward strand, genomic
                { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },             // forward strand, native
                { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 },    // forward strand, genomic, aligned
                { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 },    // forward strand, native, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },                            // forward strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },                            // forward strand, native, aligned, clipped
                { 50, 50, 50, 30, 10, 20, 10, 20, 10, 20, 10, 40, 40, 40 },             // reverse strand, genomic
                { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },             // reverse strand, native
                { 50, 50, 50, 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10, 40, 40, 40 },    // reverse strand, genomic, aligned
                { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 },    // reverse strand, native, aligned
                { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 },                            // reverse strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 }                             // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H4=3D4=3H");
        tests::CheckFrameTagsClippedAndAligned(
            "2H4=3D4=3H",                                       // CIGAR
            { 10, 20, 10, 20, 10, 20, 10, 30 },                 // input
            {
                { 10, 20, 10, 20, 10, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 10, 20, 10, 20, 10, 30 },             // forward strand, native
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 10, 20, 10, 20, 10 },             // reverse strand, genomic
                { 10, 20, 10, 20, 10, 20, 10, 30 },             // reverse strand, native
                { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },    // reverse strand, native, aligned
                { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 },    // reverse strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 }     // reverse strand, native, aligned, clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 2H3S4=3D4=3S3H");
        tests::CheckFrameTagsClippedAndAligned(
            "2H3S4=3D4=3S3H",                                                           // CIGAR
            { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },                 // input
            {
                { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },             // forward strand, genomic
                { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },             // forward strand, native
                { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 },    // forward strand, genomic, aligned
                { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 },    // forward strand, native, aligned
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },                            // forward strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 },                            // forward strand, native, aligned, clipped
                { 50, 50, 50, 30, 10, 20, 10, 20, 10, 20, 10, 40, 40, 40 },             // reverse strand, genomic
                { 40, 40, 40, 10, 20, 10, 20, 10, 20, 10, 30, 50, 50, 50 },             // reverse strand, native
                { 50, 50, 50, 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10, 40, 40, 40 },    // reverse strand, genomic, aligned
                { 40, 40, 40, 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30, 50, 50, 50 },    // reverse strand, native, aligned
                { 30, 10, 20, 10, 0, 0, 0, 20, 10, 20, 10 },                            // reverse strand, genomic, aligned, clipped
                { 10, 20, 10, 20, 0, 0, 0, 10, 20, 10, 30 }                             // reverse strand, native, aligned, clipped
            }
        );
    }
}

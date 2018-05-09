// Author: Derek Barnett

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/BamRecord.h>
#include <pbbam/BamTagCodec.h>
#include "../src/MemoryUtils.h"

// clang-format off

using namespace PacBio;
using namespace PacBio::BAM;

namespace BamRecordTests {

static
BamRecordImpl CreateBamImpl()
{
    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};

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
BamRecord CreateBam()
{ return BamRecord{ CreateBamImpl() }; }

static
void CheckRawData(const BamRecordImpl& bam)
{
    // ensure raw data (lengths at least) matches API-facing data
    const uint32_t expectedNameBytes = bam.Name().size() + 1;  // include NULL term
    const uint32_t expectedNameNulls = 4 - (expectedNameBytes % 4);
    const uint32_t expectedNameLength = expectedNameBytes + expectedNameNulls;
    const uint32_t expectedNumCigarOps = bam.CigarData().size();
    const int32_t expectedSeqLength = bam.Sequence().length();
    const size_t expectedTagsLength = BamTagCodec::Encode(bam.Tags()).size();

    //  Name        CIGAR         Sequence       Quals      Tags
    // l_qname + (n_cigar * 4) + (l_qseq+1)/2 + l_qseq + <encoded length>
    const int expectedTotalDataLength = expectedNameLength + (expectedNumCigarOps * 4) +
                                        (expectedSeqLength + 1) / 2 + expectedSeqLength +
                                        expectedTagsLength;

    const auto rawData = PacBio::BAM::internal::BamRecordMemory::GetRawData(bam);
    ASSERT_TRUE(static_cast<bool>(rawData));

    EXPECT_EQ(expectedNameNulls, rawData->core.l_extranul);
    EXPECT_EQ(expectedNameLength, rawData->core.l_qname);
    EXPECT_EQ(expectedNumCigarOps, rawData->core.n_cigar);
    EXPECT_EQ(expectedSeqLength, rawData->core.l_qseq);
    EXPECT_EQ(expectedTotalDataLength, rawData->l_data);
}

static inline
void CheckRawData(const BamRecord& bam)
{ CheckRawData(bam.Impl()); }

static
BamRecordImpl MakeCigaredImpl(const std::string& seq,
                              const std::string& cigar,
                              const Strand strand)
{
    BamRecordImpl impl;
    impl.SetMapped(true).ReferenceId(0).Position(0).MapQuality(0);
    impl.CigarData(Cigar::FromStdString(cigar));
    impl.MateReferenceId(-1).MatePosition(-1).InsertSize(0);
    impl.SetSequenceAndQualities(seq, std::string(seq.size(), '*'));
    impl.SetReverseStrand(strand == Strand::REVERSE);
    return impl;
}

static inline
BamRecord MakeCigaredRecord(const std::string& seq,
                            const std::string& cigar,
                            const Strand strand)
{ return BamRecord{ MakeCigaredImpl(seq, cigar, strand) }; }

static
BamRecord MakeCigaredBaseRecord(const std::string& bases,
                                const std::string& cigar,
                                const Strand strand)
{
    TagCollection tags;
    tags["dt"] = bases;
    tags["st"] = bases;

    const std::string seq = std::string(bases.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredFrameRecord(const std::vector<uint16_t>& frames,
                                 const std::string& cigar,
                                 const Strand strand)
{
    TagCollection tags;
    tags["ip"] = frames;
    tags["pw"] = frames;

    const std::string seq = std::string(frames.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredQualRecord(const std::string& quals,
                                const std::string& cigar,
                                const Strand strand)
{
    TagCollection tags;
    tags["dq"] = quals;
    tags["iq"] = quals;
    tags["mq"] = quals;
    tags["sq"] = quals;

    const std::string seq = std::string(quals.size(), 'N');
    BamRecordImpl impl = MakeCigaredImpl(seq, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredPulseBaseRecord(const std::string& seqBases,
                                     const std::string& pulseCalls,
                                     const std::string& pulseBases,
                                     const std::string& cigar,
                                     const Strand strand)
{
    TagCollection tags;
    tags["pc"] = pulseCalls; // PulseCall
    tags["pt"] = pulseBases; // AltLabelTag

    BamRecordImpl impl = MakeCigaredImpl(seqBases, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredPulseQualRecord(const std::string& seqBases,
                                     const std::string& pulseCalls,
                                     const std::string& pulseQuals,
                                     const std::string& cigar,
                                     const Strand strand)
{
    TagCollection tags;
    tags["pc"] = pulseCalls;
    tags["pv"] = pulseQuals; // AltLabelQV
    tags["pq"] = pulseQuals; // LabelQV
    tags["pg"] = pulseQuals; // PulseMergeQV

    BamRecordImpl impl = MakeCigaredImpl(seqBases, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredPulseFrameRecord(const std::string& seqBases,
                                     const std::string& pulseCalls,
                                     const std::vector<uint16_t>& pulseFrames,
                                     const std::string& cigar,
                                     const Strand strand)
{
    TagCollection tags;
    tags["pc"] = pulseCalls;
    tags["pd"] = pulseFrames; // PrePulseFrames
    tags["px"] = pulseFrames; // PulseCallWidth

    BamRecordImpl impl = MakeCigaredImpl(seqBases, cigar, strand);
    impl.Tags(tags);
    return BamRecord(std::move(impl));
}

static
BamRecord MakeCigaredPulseUIntRecord(const std::string& seqBases,
                                     const std::string& pulseCalls,
                                     const std::vector<uint32_t>& pulseUInts,
                                     const std::string& cigar,
                                     const Strand strand)
{
    TagCollection tags;
    tags["pc"] = pulseCalls;
    tags["sf"] = pulseUInts; // StartFrame

    BamRecordImpl impl = MakeCigaredImpl(seqBases, cigar, strand);
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

    T ForwardGenomic() const               { return d_.at(0); }
    T ForwardNative() const                { return d_.at(1); }
    T ForwardGenomicAligned() const        { return d_.at(2); }
    T ForwardNativeAligned() const         { return d_.at(3); }
    T ForwardGenomicAlignedClipped() const { return d_.at(4); }
    T ForwardNativeAlignedClipped() const  { return d_.at(5); }
    T ReverseGenomic() const               { return d_.at(6); }
    T ReverseNative() const                { return d_.at(7); }
    T ReverseGenomicAligned() const        { return d_.at(8); }
    T ReverseNativeAligned() const         { return d_.at(9); }
    T ReverseGenomicAlignedClipped() const { return d_.at(10); }
    T ReverseNativeAlignedClipped() const  { return d_.at(11); }

private:
    std::vector<T> d_;
};

// generic data type checker on the various requested states
template<typename DataType, typename MakeRecordType, typename FetchDataType>
void CheckAlignAndClip(const std::string& cigar,
                       const DataType& input,
                       const BamRecordTests::ExpectedResult<DataType>& e,
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

template<typename DataType, typename MakeRecordType, typename FetchDataType>
void CheckPulseDataAlignAndClip(const std::string& cigar,
                                const std::string& seqBases,
                                const std::string& pulseCalls,
                                const DataType& input,
                                const BamRecordTests::ExpectedResult<DataType>& allPulses,
                                const BamRecordTests::ExpectedResult<DataType>& basecallsOnly,
                                const MakeRecordType& makeRecord,
                                const FetchDataType& fetchData)
{
    {   // map to forward strand
        const BamRecord b = makeRecord(seqBases, pulseCalls, input, cigar, Strand::FORWARD);

        EXPECT_EQ(allPulses.ForwardGenomic(),               fetchData(b, Orientation::GENOMIC, false, false, PulseBehavior::ALL));
        EXPECT_EQ(allPulses.ForwardNative(),                fetchData(b, Orientation::NATIVE,  false, false, PulseBehavior::ALL));
        // no align/clipping operations available on ALL pulses

        EXPECT_EQ(basecallsOnly.ForwardGenomic(),               fetchData(b, Orientation::GENOMIC, false, false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ForwardNative(),                fetchData(b, Orientation::NATIVE,  false, false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ForwardGenomicAligned(),        fetchData(b, Orientation::GENOMIC, true,  false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ForwardNativeAligned(),         fetchData(b, Orientation::NATIVE,  true,  false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ForwardGenomicAlignedClipped(), fetchData(b, Orientation::GENOMIC, true,  true,  PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ForwardNativeAlignedClipped(),  fetchData(b, Orientation::NATIVE,  true,  true,  PulseBehavior::BASECALLS_ONLY));
    }
    {   // map to reverse strand
        const BamRecord b = makeRecord(seqBases, pulseCalls, input, cigar, Strand::REVERSE);

        EXPECT_EQ(allPulses.ReverseGenomic(),               fetchData(b, Orientation::GENOMIC, false, false, PulseBehavior::ALL));
        EXPECT_EQ(allPulses.ReverseNative(),                fetchData(b, Orientation::NATIVE,  false, false, PulseBehavior::ALL));
        // no align/clipping operations available on ALL pulses

        EXPECT_EQ(basecallsOnly.ReverseGenomic(),               fetchData(b, Orientation::GENOMIC, false, false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ReverseNative(),                fetchData(b, Orientation::NATIVE,  false, false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ReverseGenomicAligned(),        fetchData(b, Orientation::GENOMIC, true,  false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ReverseNativeAligned(),         fetchData(b, Orientation::NATIVE,  true,  false, PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ReverseGenomicAlignedClipped(), fetchData(b, Orientation::GENOMIC, true,  true,  PulseBehavior::BASECALLS_ONLY));
        EXPECT_EQ(basecallsOnly.ReverseNativeAlignedClipped(),  fetchData(b, Orientation::NATIVE,  true,  true,  PulseBehavior::BASECALLS_ONLY));
    }
}

static
void CheckBaseTagsClippedAndAligned(const std::string& cigar,
                                    const std::string& input,
                                    const ExpectedResult<std::string>& e)
{
    // aligned record + DeletionTag, SubstitutionTag
    auto makeRecord = [](const std::string& newBases,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return MakeCigaredBaseRecord(newBases, newCigar, newStrand); };

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
void CheckFrameTagsClippedAndAligned(const std::string& cigar,
                                     const std::vector<uint16_t>& input,
                                     const ExpectedResult<std::vector<uint16_t> >& e)
{

    // aligned record + IPD, PulseWidth
    auto makeRecord = [](const std::vector<uint16_t>& newFrames,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return BamRecordTests::MakeCigaredFrameRecord(newFrames, newCigar, newStrand); };

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
void CheckQualityTagsClippedAndAligned(const std::string& cigar,
                                       const std::string& input,
                                       const ExpectedResult<std::string>& e)
{
    // aligned record + DeletionQV, InsertionQV, MergeQV, SubstitutionQV
    auto makeRecord = [](const std::string& newQuals,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return BamRecordTests::MakeCigaredQualRecord(newQuals, newCigar, newStrand); };

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
void CheckQualitiesClippedAndAligned(const std::string& cigar,
                                     const std::string& input,
                                     const ExpectedResult<std::string>& e)
{
    // aligned record w/ dummy SEQ & QUALs under test
    auto makeRecord = [](const std::string& newQuals,
                         const std::string& newCigar,
                         const Strand newStrand)
    {
        const std::string seq = std::string(newQuals.size(), 'N');
        auto record = BamRecordTests::MakeCigaredRecord(seq, newCigar, newStrand);
        record.Impl().SetSequenceAndQualities(seq, newQuals);
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
void CheckSequenceClippedAndAligned(const std::string& cigar,
                                    const std::string& input,
                                    const ExpectedResult<std::string>& e)
{
    // aligned record w/ SEQ
    auto makeRecord = [](const std::string& newSeq,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return BamRecordTests::MakeCigaredRecord(newSeq, newCigar, newStrand); };

    // SEQ
    CheckAlignAndClip(cigar, input, e, makeRecord,
                      [](const BamRecord& b,
                         Orientation orientation,
                         bool aligned,
                         bool exciseSoftClips)
                      { return b.Sequence(orientation, aligned, exciseSoftClips); }
    );
}

static
void CheckPulseBaseTags(const std::string& cigar,
                        const std::string& seqBases,
                        const std::string& pulseCalls,
                        const std::string& pulseBases,
                        const ExpectedResult<std::string>& allPulses,
                        const ExpectedResult<std::string>& basecallsOnly)
{
    // aligned record + AltLabelTag
    auto makeRecord = [](const std::string& newSeqBases,
                         const std::string& newPulseCalls,
                         const std::string& newPulseBases,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return MakeCigaredPulseBaseRecord(newSeqBases, newPulseCalls, newPulseBases, newCigar, newStrand); };

    // AltLabelTag
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseBases, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.AltLabelTag(orientation, aligned, exciseSoftClips, pulseBehavior); }
    );
    // PulseCall
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseBases, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.PulseCall(orientation, aligned, exciseSoftClips, pulseBehavior); }
    );
}

static
void CheckPulseFrameTags(const std::string& cigar,
                         const std::string& seqBases,
                         const std::string& pulseCalls,
                         const std::vector<uint16_t>& pulseFrames,
                         const ExpectedResult<std::vector<uint16_t>>& allPulses,
                         const ExpectedResult<std::vector<uint16_t>>& basecallsOnly)
{
    // aligned record + PrePulseFrames
    auto makeRecord = [](const std::string& newSeqBases,
                         const std::string& newPulseCalls,
                         const std::vector<uint16_t>& newPulseFrames,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return MakeCigaredPulseFrameRecord(newSeqBases, newPulseCalls, newPulseFrames, newCigar, newStrand); };

    // PrePulseFrame
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseFrames, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.PrePulseFrames(orientation, aligned, exciseSoftClips, pulseBehavior).Data(); }
    );
    // PulseCallWidth
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseFrames, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.PulseCallWidth(orientation, aligned, exciseSoftClips, pulseBehavior).Data(); }
    );
}

/*

    { BamRecordTag::PKMEAN,            {"pa", true}  },   photons (vector<float>
    { BamRecordTag::PKMEAN_2,          {"ps", true}  },   photons
    { BamRecordTag::PKMID,             {"pm", true}  },   photons
    { BamRecordTag::PKMID_2,           {"pi", true}  },   photons
*/

static
void CheckPulseQualityTags(const std::string& cigar,
                           const std::string& seqBases,
                           const std::string& pulseCalls,
                           const std::string& pulseQuals,
                           const ExpectedResult<std::string>& allPulses,
                           const ExpectedResult<std::string>& basecallsOnly)
{
    // aligned record + AltLabelQV
    auto makeRecord = [](const std::string& newSeqBases,
                         const std::string& newPulseCalls,
                         const std::string& newPulseQuals,
                         const std::string& newCigar,
                         const Strand newStrand)
    { return MakeCigaredPulseQualRecord(newSeqBases, newPulseCalls, newPulseQuals, newCigar, newStrand); };

    // AltLabelQV
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseQuals, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.AltLabelQV(orientation, aligned, exciseSoftClips, pulseBehavior).Fastq(); }
    );
    // LabelQV
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseQuals, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.LabelQV(orientation, aligned, exciseSoftClips, pulseBehavior).Fastq(); }
    );
    // PulseMergeQV
    CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, pulseQuals, allPulses, basecallsOnly, makeRecord,
                              [](const BamRecord& b,
                                 Orientation orientation,
                                 bool aligned,
                                 bool exciseSoftClips,
                                 PulseBehavior pulseBehavior)
                              { return b.PulseMergeQV(orientation, aligned, exciseSoftClips, pulseBehavior).Fastq(); }
    );
}

static
void CheckPulseUIntTags(const std::string& cigar,
                        const std::string& seqBases,
                        const std::string& pulseCalls,
                        const std::vector<uint32_t>& startFrames,
                        const ExpectedResult<std::vector<uint32_t>>& allPulses,
                        const ExpectedResult<std::vector<uint32_t>>& basecallsOnly)
{
   // aligned record + StartFrame
   auto makeRecord = [](const std::string& newSeqBases,
                        const std::string& newPulseCalls,
                        const std::vector<uint32_t>& newStartFrames,
                        const std::string& newCigar,
                        const Strand newStrand)
   { return MakeCigaredPulseUIntRecord(newSeqBases, newPulseCalls, newStartFrames, newCigar, newStrand); };

   // StartFrame
   CheckPulseDataAlignAndClip(cigar, seqBases, pulseCalls, startFrames, allPulses, basecallsOnly, makeRecord,
                             [](const BamRecord& b,
                                Orientation orientation,
                                bool aligned,
                                bool exciseSoftClips,
                                PulseBehavior pulseBehavior)
                             { return b.StartFrame(orientation, aligned, exciseSoftClips, pulseBehavior); }
   );
}



} // namespace BamRecordTests

TEST(BamRecordTest, DefaultValues)
{
    BamRecord bam;
    const std::string emptyString;

    // BamRecordImpl data
    EXPECT_EQ(0, bam.Impl().Bin());
    EXPECT_EQ(BamRecordImpl::UNMAPPED, bam.Impl().Flag());  // forced init unmapped
    EXPECT_EQ(0, bam.Impl().InsertSize());
    EXPECT_EQ(255, bam.Impl().MapQuality());
    EXPECT_EQ(-1, bam.Impl().MateReferenceId());
    EXPECT_EQ(-1, bam.Impl().MatePosition());
    EXPECT_EQ(-1, bam.Impl().Position());
    EXPECT_EQ(-1, bam.Impl().ReferenceId());
    EXPECT_EQ(0, bam.Impl().Tags().size());

    EXPECT_FALSE(bam.Impl().IsDuplicate());
    EXPECT_FALSE(bam.Impl().IsFailedQC());
    EXPECT_FALSE(bam.Impl().IsFirstMate());
    EXPECT_FALSE(bam.Impl().IsMapped());             // forced init unmapped
    EXPECT_TRUE(bam.Impl().IsMateMapped());
    EXPECT_FALSE(bam.Impl().IsMateReverseStrand());
    EXPECT_FALSE(bam.Impl().IsPaired());
    EXPECT_TRUE(bam.Impl().IsPrimaryAlignment());
    EXPECT_FALSE(bam.Impl().IsProperPair());
    EXPECT_FALSE(bam.Impl().IsReverseStrand());
    EXPECT_FALSE(bam.Impl().IsSecondMate());
    EXPECT_FALSE(bam.Impl().IsSupplementaryAlignment());

    EXPECT_EQ(emptyString, bam.Impl().Name());
    EXPECT_EQ(emptyString, bam.Impl().CigarData().ToStdString());
    EXPECT_EQ(emptyString, bam.Impl().Sequence());
    EXPECT_EQ(emptyString, bam.Impl().Qualities().Fastq());

    // PacBio data
    EXPECT_EQ(-1, bam.AlignedStart());
    EXPECT_EQ(-1, bam.AlignedEnd());

    EXPECT_FALSE(bam.HasHoleNumber());
    EXPECT_FALSE(bam.HasNumPasses());
    EXPECT_FALSE(bam.HasQueryEnd());
    EXPECT_FALSE(bam.HasQueryStart());
    EXPECT_FALSE(bam.HasReadAccuracy());

    EXPECT_THROW(bam.HoleNumber(), std::exception);
    EXPECT_THROW(bam.NumPasses(), std::exception);
    EXPECT_EQ(int32_t{0}, bam.QueryEnd());
    EXPECT_EQ(int32_t{0}, bam.QueryStart());
    EXPECT_THROW(bam.ReadAccuracy(), std::exception);

    EXPECT_FALSE(bam.HasDeletionQV());
    EXPECT_FALSE(bam.HasDeletionTag());
    EXPECT_FALSE(bam.HasInsertionQV());
    EXPECT_FALSE(bam.HasMergeQV());
    EXPECT_FALSE(bam.HasSubstitutionQV());
    EXPECT_FALSE(bam.HasSubstitutionTag());

    EXPECT_THROW(bam.DeletionQV(),      std::exception);
    EXPECT_THROW(bam.DeletionTag(),     std::exception);
    EXPECT_THROW(bam.InsertionQV(),     std::exception);
    EXPECT_THROW(bam.MergeQV(),         std::exception);
    EXPECT_THROW(bam.SubstitutionQV(),  std::exception);
    EXPECT_THROW(bam.SubstitutionTag(), std::exception);

    // raw data
    BamRecordTests::CheckRawData(bam);
}

TEST(BamRecordTest, FromBamRecordImpl)
{
    // check generic data
    BamRecordImpl genericBam = BamRecordTests::CreateBamImpl();

    EXPECT_EQ(42, genericBam.Bin());
    EXPECT_EQ(42, genericBam.Flag());
    EXPECT_EQ(42, genericBam.InsertSize());
    EXPECT_EQ(42, genericBam.MapQuality());
    EXPECT_EQ(42, genericBam.MateReferenceId());
    EXPECT_EQ(42, genericBam.MatePosition());
    EXPECT_EQ(42, genericBam.Position());
    EXPECT_EQ(42, genericBam.ReferenceId());

    const TagCollection genericTags = genericBam.Tags();
    EXPECT_TRUE(genericTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), genericTags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, genericTags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), genericTags.at("CA").ToUInt8Array());

    // copy ctor
    BamRecord bam1(genericBam);

    EXPECT_EQ(42, bam1.Impl().Bin());
    EXPECT_EQ(42, bam1.Impl().Flag());
    EXPECT_EQ(42, bam1.Impl().InsertSize());
    EXPECT_EQ(42, bam1.Impl().MapQuality());
    EXPECT_EQ(42, bam1.Impl().MateReferenceId());
    EXPECT_EQ(42, bam1.Impl().MatePosition());
    EXPECT_EQ(42, bam1.Impl().Position());
    EXPECT_EQ(42, bam1.Impl().ReferenceId());

    const TagCollection bam1Tags = bam1.Impl().Tags();
    EXPECT_TRUE(bam1Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), bam1Tags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, bam1Tags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), bam1Tags.at("CA").ToUInt8Array());

    // copy assignment
    BamRecord bam2;
    bam2 = genericBam;

    EXPECT_EQ(42, bam2.Impl().Bin());
    EXPECT_EQ(42, bam2.Impl().Flag());
    EXPECT_EQ(42, bam2.Impl().InsertSize());
    EXPECT_EQ(42, bam2.Impl().MapQuality());
    EXPECT_EQ(42, bam2.Impl().MateReferenceId());
    EXPECT_EQ(42, bam2.Impl().MatePosition());
    EXPECT_EQ(42, bam2.Impl().Position());
    EXPECT_EQ(42, bam2.Impl().ReferenceId());

    const TagCollection bam2Tags = bam2.Impl().Tags();
    EXPECT_TRUE(bam2Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), bam2Tags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, bam2Tags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), bam2Tags.at("CA").ToUInt8Array());

    // change genericBam, make sure we deep copied bam1 & bam2
    genericBam.Position(2000);

    EXPECT_EQ(2000, genericBam.Position());
    EXPECT_EQ(42, bam1.Impl().Position());
    EXPECT_EQ(42, bam2.Impl().Position());

    // move ctor
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    BamRecord bam3(std::move(BamRecordTests::CreateBamImpl()));
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    EXPECT_EQ(42, bam3.Impl().Bin());
    EXPECT_EQ(42, bam3.Impl().Flag());
    EXPECT_EQ(42, bam3.Impl().InsertSize());
    EXPECT_EQ(42, bam3.Impl().MapQuality());
    EXPECT_EQ(42, bam3.Impl().MateReferenceId());
    EXPECT_EQ(42, bam3.Impl().MatePosition());
    EXPECT_EQ(42, bam3.Impl().Position());
    EXPECT_EQ(42, bam3.Impl().ReferenceId());

    const TagCollection bam3Tags = bam3.Impl().Tags();
    EXPECT_TRUE(bam3Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), bam3Tags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, bam3Tags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), bam3Tags.at("CA").ToUInt8Array());

    // move assignment
    BamRecord bam4;
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
    bam4 = std::move(BamRecordTests::CreateBamImpl());
#ifdef __clang__
#pragma clang diagnostic pop
#endif

    EXPECT_EQ(42, bam4.Impl().Bin());
    EXPECT_EQ(42, bam4.Impl().Flag());
    EXPECT_EQ(42, bam4.Impl().InsertSize());
    EXPECT_EQ(42, bam4.Impl().MapQuality());
    EXPECT_EQ(42, bam4.Impl().MateReferenceId());
    EXPECT_EQ(42, bam4.Impl().MatePosition());
    EXPECT_EQ(42, bam4.Impl().Position());
    EXPECT_EQ(42, bam4.Impl().ReferenceId());

    const TagCollection bam4Tags = bam4.Impl().Tags();
    EXPECT_TRUE(bam4Tags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), bam4Tags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, bam4Tags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), bam4Tags.at("CA").ToUInt8Array());
}

TEST(BamRecordTest, SelfAssignmentTolerated)
{
    BamRecord bam1;
    bam1.Impl().Bin(42);
    bam1.Impl().Flag(42);
    bam1.Impl().InsertSize(42);
    bam1.Impl().MapQuality(42);
    bam1.Impl().MatePosition(42);
    bam1.Impl().MateReferenceId(42);
    bam1.Impl().Position(42);
    bam1.Impl().ReferenceId(42);

    TagCollection tags;
    tags["HX"] = std::string("1abc75");
    tags["HX"].Modifier(TagModifier::HEX_STRING);
    tags["CA"] = std::vector<uint8_t>({34, 5, 125});
    tags["XY"] = int32_t{-42};
    bam1.Impl().Tags(tags);

    bam1 = bam1;

    EXPECT_EQ(42, bam1.Impl().Bin());
    EXPECT_EQ(42, bam1.Impl().Flag());
    EXPECT_EQ(42, bam1.Impl().InsertSize());
    EXPECT_EQ(42, bam1.Impl().MapQuality());
    EXPECT_EQ(42, bam1.Impl().MateReferenceId());
    EXPECT_EQ(42, bam1.Impl().MatePosition());
    EXPECT_EQ(42, bam1.Impl().Position());
    EXPECT_EQ(42, bam1.Impl().ReferenceId());

    const TagCollection fetchedTags1 = bam1.Impl().Tags();
    EXPECT_TRUE(fetchedTags1.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags1.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags1.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags1.at("CA").ToUInt8Array());

    BamRecordTests::CheckRawData(bam1);
}

TEST(BamRecordTest, CoreSetters)
{
    // create basic BAM with (generic) data
    BamRecord bam = BamRecordTests::CreateBam();

    QualityValues testQVs;
    testQVs.push_back(0);
    testQVs.push_back(1);

    const std::string testTags = "GATTACA";

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
    EXPECT_EQ(42, bam.Impl().Bin());
    EXPECT_EQ(42, bam.Impl().Flag());
    EXPECT_EQ(42, bam.Impl().InsertSize());
    EXPECT_EQ(42, bam.Impl().MapQuality());
    EXPECT_EQ(42, bam.Impl().MateReferenceId());
    EXPECT_EQ(42, bam.Impl().MatePosition());
    EXPECT_EQ(42, bam.Impl().Position());
    EXPECT_EQ(42, bam.Impl().ReferenceId());

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
    const TagCollection fetchedTags = bam.Impl().Tags();
    EXPECT_TRUE(fetchedTags.at("HX").HasModifier(TagModifier::HEX_STRING));
    EXPECT_EQ(std::string("1abc75"), fetchedTags.at("HX").ToString());
    EXPECT_EQ(int32_t{-42}, fetchedTags.at("XY").ToInt32());
    EXPECT_EQ(std::vector<uint8_t>({34, 5, 125}), fetchedTags.at("CA").ToUInt8Array());

    BamRecordTests::CheckRawData(bam);
}

TEST(BamRecordTest, SequenceOrientation)
{
    {
        SCOPED_TRACE("Simple CIGAR Sequence");
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckQualitiesClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
            "3=4N3=",       // CIGAR
            "ACGTTT",        // input
            {
                "ACGTTT",    // forward strand, genomic
                "ACGTTT",    // forward strand, native
                "ACGTTT",    // forward strand, genomic, aligned
                "ACGTTT",    // forward strand, native,  aligned
                "ACGTTT",    // forward strand, genomic, aligned + clipped
                "ACGTTT",    // forward strand, native,  aligned + clipped
                "ACGTTT",    // reverse strand, genomic
                "AAACGT",    // reverse strand, native
                "ACGTTT",    // reverse strand, genomic, aligned
                "AAACGT",    // reverse strand, native,  aligned
                "ACGTTT",    // reverse strand, genomic, aligned + clipped
                "AAACGT"     // reverse strand, native,  aligned + clipped
            }
        );
    }
    {
        SCOPED_TRACE("CIGAR: 1S8=1S");
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckSequenceClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckQualityTagsClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckBaseTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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
        BamRecordTests::CheckFrameTagsClippedAndAligned(
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

TEST(BamRecordTest, PulseBaseTags)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        BamRecordTests::CheckPulseBaseTags(
            "4=3D4=",           // CIGAR
            "AACCGTTA",         // seqBases
            "AAaaCCGggTTA",     // pulseCalls
            "AAaaCCGggTTA",     // tag data

            {   // all pulses

                "AAaaCCGggTTA",     // forward strand, genomic
                "AAaaCCGggTTA",     // forward strand, native
                "",  // forward strand, genomic, aligned
                "",  // forward strand, native, aligned
                "",  // forward strand, genomic, aligned, clipped
                "",  // forward strand, native, aligned, clipped
                "TAAccCGGttTT",     // reverse strand, genomic
                "AAaaCCGggTTA",     // reverse strand, native
                "",  // reverse strand, genomic, aligned
                "",  // reverse strand, native, aligned
                "",  // reverse strand, genomic, aligned, clipped
                ""   // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseBaseTags(
            "4=1D2I2D4=",       // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            "ATttCCTtAGGggTT",  // tag data

            {   // all pulses

                "ATttCCTtAGGggTT",       // forward strand, genomic
                "ATttCCTtAGGggTT",       // forward strand, native
                "",    // forward strand, genomic, aligned
                "",    // forward strand, native, aligned
                "",    // forward strand, genomic, aligned, clipped
                "",    // forward strand, native, aligned, clipped
                "AAccCCTaAGGaaAT",       // reverse strand, genomic
                "ATttCCTtAGGggTT",       // reverse strand, native
                "",    // reverse strand, genomic, aligned
                "",    // reverse strand, native, aligned
                "",    // reverse strand, genomic, aligned, clipped
                ""     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseBaseTags(
            "4=1D2P2I2P2D4=",   // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            "ATttCCTtAGGggTT",  // tag data
            {
                "ATttCCTtAGGggTT",           // forward strand, genomic
                "ATttCCTtAGGggTT",           // forward strand, native
                "",    // forward strand, genomic, aligned
                "",    // forward strand, native, aligned
                "",    // forward strand, genomic, aligned, clipped
                "",    // forward strand, native, aligned, clipped
                "AAccCCTaAGGaaAT",           // reverse strand, genomic
                "ATttCCTtAGGggTT",           // reverse strand, native
                "",    // reverse strand, genomic, aligned
                "",    // reverse strand, native, aligned
                "",    // reverse strand, genomic, aligned, clipped
                ""     // reverse strand, native, aligned, clipped
            },
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
        BamRecordTests::CheckPulseBaseTags(
            "3S4=3D4=3S",               // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            "TTTttAACCccGTTAaaCCG",     // tag data

            {   // all pulses

                "TTTttAACCccGTTAaaCCG",       // forward strand, genomic
                "TTTttAACCccGTTAaaCCG",       // forward strand, native
                "",         // forward strand, genomic, aligned
                "",         // forward strand, native, aligned
                "",          // forward strand, genomic, aligned, clipped
                "",          // forward strand, native, aligned, clipped
                "CGGttTAACggGGTTaaAAA",       // reverse strand, genomic
                "TTTttAACCccGTTAaaCCG",       // reverse strand, native
                "",    // reverse strand, genomic, aligned
                "",    // reverse strand, native, aligned
                "",     // reverse strand, genomic, aligned, clipped
                ""           // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseBaseTags(
            "2H4=3D4=3H",       // CIGAR
            "AACCGTTA",         // seqBases
            "AAaaCCGggTTA",     // pulseCalls
            "AAaaCCGggTTA",     // tag data

            {   // all pulses

                "AAaaCCGggTTA",     // forward strand, genomic
                "AAaaCCGggTTA",     // forward strand, native
                "",  // forward strand, genomic, aligned
                "",  // forward strand, native, aligned
                "",  // forward strand, genomic, aligned, clipped
                "",  // forward strand, native, aligned, clipped
                "TAAccCGGttTT",     // reverse strand, genomic
                "AAaaCCGggTTA",     // reverse strand, native
                "",  // reverse strand, genomic, aligned
                "",  // reverse strand, native, aligned
                "",  // reverse strand, genomic, aligned, clipped
                ""   // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseBaseTags(
            "2H3S4=3D4=3S3H",           // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            "TTTttAACCccGTTAaaCCG",     // tag data

            {   // all pulses

                "TTTttAACCccGTTAaaCCG",       // forward strand, genomic
                "TTTttAACCccGTTAaaCCG",       // forward strand, native
                "",         // forward strand, genomic, aligned
                "",         // forward strand, native, aligned
                "",          // forward strand, genomic, aligned, clipped
                "",          // forward strand, native, aligned, clipped
                "CGGttTAACggGGTTaaAAA",       // reverse strand, genomic
                "TTTttAACCccGTTAaaCCG",       // reverse strand, native
                "",         // reverse strand, genomic, aligned
                "",         // reverse strand, native, aligned
                "",          // reverse strand, genomic, aligned, clipped
                ""           // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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

TEST(BamRecordTest, PulseQualityTags)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        BamRecordTests::CheckPulseQualityTags(
            "4=3D4=",           // CIGAR
            "AACCGTTA",         // seqBases
            "AAaaCCGggTTA",     // pulseCalls
            "?]!!?]?!!]?@",     // tag data

            {   // all pulses

                "?]!!?]?!!]?@",     // forward strand, genomic
                "?]!!?]?!!]?@",     // forward strand, native
                "",  // forward strand, genomic, aligned
                "",  // forward strand, native,  aligned
                "",  // forward strand, genomic, aligned + clipped
                "",  // forward strand, native,  aligned + clipped
                "@?]!!?]?!!]?",     // reverse strand, genomic
                "?]!!?]?!!]?@",     // reverse strand, native
                "",  // reverse strand, genomic, aligned
                "",  // reverse strand, native,  aligned
                "",  // reverse strand, genomic, aligned + clipped
                ""   // reverse strand, native,  aligned + clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseQualityTags(
            "4=1D2I2D4=",       // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            "?]!!?]8!7?]!!?@",  // tag data

            {   // all pulses

                "?]!!?]8!7?]!!?@",       // forward strand, genomic
                "?]!!?]8!7?]!!?@",       // forward strand, native
                "",    // forward strand, genomic, aligned
                "",    // forward strand, native,  aligned
                "",    // forward strand, genomic, aligned + clipped
                "",    // forward strand, native,  aligned + clipped
                "@?!!]?7!8]?!!]?",       // reverse strand, genomic
                "?]!!?]8!7?]!!?@",       // reverse strand, native
                "",    // reverse strand, genomic, aligned
                "",    // reverse strand, native,  aligned
                "",    // reverse strand, genomic, aligned + clipped
                ""     // reverse strand, native,  aligned + clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseQualityTags(
            "4=1D2P2I2P2D4=",   // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            "?]!!?]8!7?]!!?@",  // tag data
        {
            "?]!!?]8!7?]!!?@",           // forward strand, genomic
            "?]!!?]8!7?]!!?@",           // forward strand, native
            "",    // forward strand, genomic, aligned
            "",    // forward strand, native,  aligned
            "",    // forward strand, genomic, aligned + clipped
            "",    // forward strand, native,  aligned + clipped
            "@?!!]?7!8]?!!]?",           // reverse strand, genomic
            "?]!!?]8!7?]!!?@",           // reverse strand, native
            "",    // reverse strand, genomic, aligned
            "",    // reverse strand, native,  aligned
            "",    // reverse strand, genomic, aligned + clipped
            ""     // reverse strand, native,  aligned + clipped
        },
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
        BamRecordTests::CheckPulseQualityTags(
            "3S4=3D4=3S",               // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            "vvv!!?]?]!!?]?@!!xxx",     // tag data

            {   // all pulses

                "vvv!!?]?]!!?]?@!!xxx",       // forward strand, genomic
                "vvv!!?]?]!!?]?@!!xxx",       // forward strand, native
                "",    // forward strand, genomic, aligned
                "",    // forward strand, native, aligned
                "",          // forward strand, genomic, aligned, clipped
                "",          // forward strand, native, aligned, clipped
                "xxx!!@?]?!!]?]?!!vvv",       // reverse strand, genomic
                "vvv!!?]?]!!?]?@!!xxx",       // reverse strand, native
                "",    // reverse strand, genomic, aligned
                "",    // reverse strand, native, aligned
                "",          // reverse strand, genomic, aligned, clipped
                ""           // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseQualityTags(
            "2H4=3D4=3H",       // CIGAR
            "AACCGTTA",         // seqBases
            "AAaaCCGggTTA",     // pulseCalls
            "?]!!?]?!!]?@",     // tag data

            {   // all pulses

                "?]!!?]?!!]?@",     // forward strand, genomic
                "?]!!?]?!!]?@",     // forward strand, native
                "",  // forward strand, genomic, aligned
                "",  // forward strand, native, aligned
                "",  // forward strand, genomic, aligned, clipped
                "",  // forward strand, native, aligned, clipped
                "@?]!!?]?!!]?",     // reverse strand, genomic
                "?]!!?]?!!]?@",     // reverse strand, native
                "",  // reverse strand, genomic, aligned
                "",  // reverse strand, native, aligned
                "",  // reverse strand, genomic, aligned, clipped
                ""   // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseQualityTags(
            "2H3S4=3D4=3S3H",           // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            "vvv!!?]?]!!?]?@!!xxx",     // tag data

            {   // all pulses

                "vvv!!?]?]!!?]?@!!xxx",       // forward strand, genomic
                "vvv!!?]?]!!?]?@!!xxx",       // forward strand, native
                "",    // forward strand, genomic, aligned
                "",    // forward strand, native, aligned
                "",          // forward strand, genomic, aligned, clipped
                "",          // forward strand, native, aligned, clipped
                "xxx!!@?]?!!]?]?!!vvv",       // reverse strand, genomic
                "vvv!!?]?]!!?]?@!!xxx",       // reverse strand, native
                "",    // reverse strand, genomic, aligned
                "",    // reverse strand, native, aligned
                "",          // reverse strand, genomic, aligned, clipped
                ""           // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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

TEST(BamRecordTest, PulseFrameTags)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        BamRecordTests::CheckPulseFrameTags(
            "4=3D4=",       // CIGAR
            "AACCGTTA",     // seqBases
            "AAaaCCGggTTA", // pulseCalls
            { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },   // tag data

            {   // all pulses

                { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseFrameTags(
            "4=1D2I2D4=",       // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 }, // tag data

            {   // all pulses

                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 0, 0, 20, 10, 70, 0, 80, 20, 10, 0, 0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseFrameTags(
            "4=1D2P2I2P2D4=",   // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 }, // tag data

            {   // all pulses

                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 0, 0, 20, 10, 70, 0, 80, 20, 10, 0, 0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseFrameTags(
            "3S4=3D4=3S",               // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },   // tag data

            {   // all pulses

                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 50, 50, 50, 0, 0, 30, 10, 20, 10, 0, 0, 20, 10, 20, 10, 0, 0, 40, 40, 40 },             // reverse strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseFrameTags(
            "2H4=3D4=3H",       // CIGAR
            "AACCGTTA",         // seqBases
            "AAaaCCGggTTA",     // pulseCalls
            { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 }, // tag data

            {   // all pulses

                { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseFrameTags(
            "2H3S4=3D4=3S3H",           // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },                 // tag data

            {   // all pulses

                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 50, 50, 50, 0, 0, 30, 10, 20, 10, 0, 0, 20, 10, 20, 10, 0, 0, 40, 40, 40 },             // reverse strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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

TEST(BamRecordTest, PulseUIntTags)
{
    {
        SCOPED_TRACE("CIGAR: 4=3D4=");
        BamRecordTests::CheckPulseUIntTags(
            "4=3D4=",       // CIGAR
            "AACCGTTA",     // seqBases
            "AAaaCCGggTTA", // pulseCalls
            { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },   // tag data

            {   // all pulses

                { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0,0, 10, 20, 10, 0,0, 20, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseUIntTags(
            "4=1D2I2D4=",       // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 }, // tag data

            {   // all pulses

                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 0, 0, 20, 10, 70, 0, 80, 20, 10, 0, 0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseUIntTags(
            "4=1D2P2I2P2D4=",   // CIGAR
            "ATCCTAGGTT",       // seqBases
            "ATttCCTtAGGggTT",  // pulseCalls
            { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 }, // tag data

            {   // all pulses

                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 0, 0, 20, 10, 70, 0, 80, 20, 10, 0, 0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0, 0, 10, 20, 80, 0, 70, 10, 20, 0, 0, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseUIntTags(
            "3S4=3D4=3S",               // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },   // tag data

            {   // all pulses

                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 50, 50, 50, 0, 0, 30, 10, 20, 10, 0, 0, 20, 10, 20, 10, 0, 0, 40, 40, 40 },             // reverse strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseUIntTags(
            "2H4=3D4=3H",       // CIGAR
            "AACCGTTA",         // seqBases
            "AAaaCCGggTTA",     // pulseCalls
            { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 }, // tag data

            {   // all pulses

                { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 },             // forward strand, genomic
                { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 30, 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10 },             // reverse strand, genomic
                { 10, 20, 0, 0, 10, 20, 10, 0, 0, 20, 10, 30 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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
        BamRecordTests::CheckPulseUIntTags(
            "2H3S4=3D4=3S3H",           // CIGAR
            "TTTAACCGTTACCG",           // seqBases
            "TTTttAACCccGTTAaaCCG",     // pulseCalls
            { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },                 // tag data

            {   // all pulses

                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // forward strand, native
                { },    // forward strand, genomic, aligned
                { },    // forward strand, native, aligned
                { },    // forward strand, genomic, aligned, clipped
                { },    // forward strand, native, aligned, clipped
                { 50, 50, 50, 0, 0, 30, 10, 20, 10, 0, 0, 20, 10, 20, 10, 0, 0, 40, 40, 40 },             // reverse strand, genomic
                { 40, 40, 40, 0, 0, 10, 20, 10, 20, 0, 0, 10, 20, 10, 30, 0, 0, 50, 50, 50 },             // reverse strand, native
                { },    // reverse strand, genomic, aligned
                { },    // reverse strand, native, aligned
                { },    // reverse strand, genomic, aligned, clipped
                { }     // reverse strand, native, aligned, clipped
            },
            {   // basecalls only

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

TEST(BamRecordTest, PulseExclusionTag)
{
    const std::vector<PacBio::BAM::PulseExclusionReason> reasons =
    {
        PulseExclusionReason::BASE
      , PulseExclusionReason::PAUSE
      , PulseExclusionReason::SHORT_PULSE
      , PulseExclusionReason::BURST
      , PulseExclusionReason::BASE
      , PulseExclusionReason::PAUSE
    };

    auto bam = BamRecordTests::CreateBam();
    bam.PulseExclusionReason(reasons);

    EXPECT_TRUE(bam.HasPulseExclusion());
    auto result = bam.PulseExclusionReason();
    EXPECT_EQ(reasons, result);

}

TEST(BamRecordTest, TranscriptRecord)
{
    const std::string readTypeStr{"TRANSCRIPT"};
    const auto readGroupId = MakeReadGroupId("transcript", readTypeStr);

    ReadGroupInfo rg{readGroupId};
    rg.ReadType(readTypeStr);

    BamHeader header;
    header.Version("1.1")
        .SortOrder("queryname")
        .PacBioBamVersion("3.0.1");

    BamRecord bam{header};
    bam.Impl().Name("transcript/1234");

    EXPECT_EQ(RecordType::TRANSCRIPT, bam.Type());
    EXPECT_EQ(1234, bam.HoleNumber());
    EXPECT_THROW({bam.QueryStart();}, std::runtime_error);
    EXPECT_THROW({bam.QueryEnd();}, std::runtime_error);
}

// clang-format on

// File Description
/// \file SimpleRead.cpp
/// \brief Implements the SimpleRead class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SimpleRead.h"

#include <cassert>

#include <stdexcept>
#include <type_traits>

#include "Clipping.h"
#include "SequenceUtils.h"
#include "SimpleReadImpl.h"

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T>
T clipContainer(const T& input, const size_t pos, const size_t len)
{
    if (input.empty()) return {};
    assert(input.size() >= pos + len);
    return T{input.cbegin() + pos, input.cbegin() + pos + len};
}

void ClipSimpleRead(SimpleRead& read, const internal::ClipResult& result, size_t start, size_t end)
{
    const auto clipFrom = result.clipOffset_;
    const auto clipLength = (end - start);
    read.Sequence = clipContainer(read.Sequence, clipFrom, clipLength);
    read.Qualities = clipContainer(read.Qualities, clipFrom, clipLength);
    read.QueryStart = result.qStart_;
    read.QueryEnd = result.qEnd_;
    if (read.PulseWidths)
        read.PulseWidths = clipContainer(read.PulseWidths->Data(), clipFrom, clipLength);
    if (read.IPD) read.IPD = clipContainer(read.IPD->Data(), clipFrom, clipLength);
}

// NOTE: 'result' is moved into here, so we can take the CIGAR
void ClipMappedRead(MappedSimpleRead& read, internal::ClipResult result)
{
    // clip common data
    ClipSimpleRead(read, result, result.qStart_, result.qEnd_);

    // clip mapped data
    read.Cigar = std::move(result.cigar_);
    read.TemplateStart = result.refPos_;
    read.TemplateEnd = read.TemplateStart + ReferenceLength(read.Cigar);
}

}  // namespace internal

//
// SimpleRead
//

static_assert(std::is_copy_constructible<SimpleRead>::value,
              "SimpleRead(const SimpleRead&) is not = default");
static_assert(std::is_copy_assignable<SimpleRead>::value,
              "SimpleRead& operator=(const SimpleRead&) is not = default");

static_assert(std::is_nothrow_move_constructible<SimpleRead>::value,
              "SimpleRead(SimpleRead&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<SimpleRead>::value ==
                  std::is_nothrow_move_assignable<std::string>::value,
              "");

SimpleRead::SimpleRead(const BamRecord& bam)
    : Name{bam.FullName()}
    , Sequence{bam.Sequence()}
    , Qualities{bam.Qualities()}
    , SignalToNoise{bam.SignalToNoise()}
    , QueryStart{bam.QueryStart()}
    , QueryEnd{bam.QueryEnd()}
{
    if (bam.IsMapped() && (bam.AlignedStrand() == Strand::REVERSE)) {
        ReverseComplement(Sequence);
        Reverse(Qualities);
    }
}

SimpleRead::SimpleRead(std::string name, std::string seq, QualityValues qualities, SNR snr)
    : Name{std::move(name)}
    , Sequence{std::move(seq)}
    , Qualities{std::move(qualities)}
    , SignalToNoise{std::move(snr)}
{
    QueryStart = 0;
    QueryEnd = Sequence.size();
}

SimpleRead::SimpleRead(std::string name, std::string seq, QualityValues qualities, SNR snr,
                       Position qStart, Position qEnd)
    : Name{std::move(name)}
    , Sequence{std::move(seq)}
    , Qualities{std::move(qualities)}
    , SignalToNoise{std::move(snr)}
    , QueryStart{qStart}
    , QueryEnd{qEnd}
{
}

SimpleRead::SimpleRead(std::string name, std::string seq, QualityValues qualities, SNR snr,
                       Position qStart, Position qEnd, Frames pulseWidths, Frames ipd)
    : Name{std::move(name)}
    , Sequence{std::move(seq)}
    , Qualities{std::move(qualities)}
    , SignalToNoise{std::move(snr)}
    , QueryStart{qStart}
    , QueryEnd{qEnd}
    , PulseWidths{std::move(pulseWidths)}
    , IPD{std::move(ipd)}
{
}

//
// MappedSimpleRead
//

static_assert(std::is_copy_constructible<MappedSimpleRead>::value,
              "MappedSimpleRead(const MappedSimpleRead&) is not = default");
static_assert(std::is_copy_assignable<MappedSimpleRead>::value,
              "MappedSimpleRead& operator=(const MappedSimpleRead&) is not = default");

static_assert(std::is_nothrow_move_constructible<MappedSimpleRead>::value,
              "MappedSimpleRead(MappedSimpleRead&&) is not = noexcept");
static_assert(std::is_nothrow_move_assignable<MappedSimpleRead>::value ==
                  std::is_nothrow_move_assignable<SimpleRead>::value,
              "");

MappedSimpleRead::MappedSimpleRead(const BamRecord& bam)
    : SimpleRead{bam}
    , Strand{bam.AlignedStrand()}
    , TemplateStart{bam.ReferenceStart()}
    , TemplateEnd{bam.ReferenceEnd()}
    , Cigar{bam.CigarData()}
    , MapQuality{bam.MapQuality()}
{
    if (!bam.IsMapped()) {
        throw std::runtime_error{"MappedSimpleRead error: input BAM record '" + bam.FullName() +
                                 "' is not mapped"};
    }
}

MappedSimpleRead::MappedSimpleRead(const SimpleRead& read, PacBio::BAM::Strand strand,
                                   Position templateStart, Position templateEnd,
                                   PacBio::BAM::Cigar cigar, uint8_t mapQV)
    : SimpleRead{read}
    , Strand{strand}
    , TemplateStart{templateStart}
    , TemplateEnd{templateEnd}
    , Cigar{std::move(cigar)}
    , MapQuality{mapQV}
{
}

//
// Clipping helpers
//

void ClipToQuery(SimpleRead& read, Position start, Position end)
{
    // skip out if clip not needed
    if (start <= read.QueryStart && end >= read.QueryEnd) return;

    // calculate clipping
    internal::ClipToQueryConfig clipConfig{read.Sequence.size(),
                                           read.QueryStart,
                                           read.QueryEnd,
                                           start,
                                           end,
                                           // TODO: these next are unused for unmapped clip-to-query
                                           UnmappedPosition,
                                           Strand::FORWARD,
                                           {},
                                           false};
    auto result = internal::ClipToQuery(clipConfig);

    // apply clipping
    internal::ClipSimpleRead(read, std::move(result), start, end);
}

void ClipToQuery(MappedSimpleRead& read, Position start, Position end)
{
    // skip out if clip not needed
    if (start <= read.QueryStart && end >= read.QueryEnd) return;

    // calculate clipping
    internal::ClipToQueryConfig clipConfig{read.Sequence.size(),
                                           read.QueryStart,
                                           read.QueryEnd,
                                           start,
                                           end,
                                           read.TemplateStart,
                                           read.Strand,
                                           std::move(read.Cigar),  // we'll get this back
                                           true};
    auto result = internal::ClipToQuery(clipConfig);

    // apply clipping
    internal::ClipMappedRead(read, std::move(result));
}

void ClipToReference(MappedSimpleRead& read, Position start, Position end,
                     bool exciseFlankingInserts)
{
    // return emptied read if clip region is disjoint from
    if (end <= read.TemplateStart || start >= read.TemplateEnd) {
        read.Sequence.clear();
        read.Qualities.clear();
        read.QueryStart = -1;
        read.QueryEnd = -1;
        if (read.PulseWidths) read.PulseWidths->DataRaw().clear();
        read.TemplateStart = -1;
        read.TemplateEnd = -1;
        read.Cigar.clear();
        read.MapQuality = 255;
        return;
    }

    // skip out if clip region covers aligned region (no clip needed)
    if (start <= read.TemplateStart && end >= read.TemplateEnd) return;

    // calculate clipping
    internal::ClipToReferenceConfig clipConfig{
        internal::ClipToQueryConfig{read.Sequence.size(), read.QueryStart, read.QueryEnd, start,
                                    end, read.TemplateStart, read.Strand,
                                    std::move(read.Cigar),  // we'll get this back
                                    true},
        read.TemplateEnd, start, end, exciseFlankingInserts};
    auto result = internal::ClipToReference(clipConfig);

    // apply clipping
    internal::ClipMappedRead(read, std::move(result));
}

}  // namespace BAM
}  // namespace PacBio

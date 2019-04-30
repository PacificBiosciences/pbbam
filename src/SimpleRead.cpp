// File Description
/// \file SimpleRead.cpp
/// \brief Implements the SimpleRead class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SimpleRead.h"

#include "Clipping.h"
#include "SequenceUtils.h"

namespace PacBio {
namespace BAM {

//
// SimpleRead
//

SimpleRead::SimpleRead(const BamRecord& bam)
    : Name{bam.FullName()}
    , Sequence{bam.Sequence()}
    , Qualities{bam.Qualities()}
    , SignalToNoise{bam.SignalToNoise()}
    , QueryStart{bam.QueryStart()}
    , QueryEnd{bam.QueryEnd()}
{
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

SimpleRead::SimpleRead(const SimpleRead&) = default;

SimpleRead::SimpleRead(SimpleRead&&) = default;

SimpleRead& SimpleRead::operator=(const SimpleRead&) = default;

SimpleRead& SimpleRead::operator=(SimpleRead&&) = default;

SimpleRead::~SimpleRead() = default;

//
// MappedSimpleRead
//

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

MappedSimpleRead::MappedSimpleRead(const MappedSimpleRead&) = default;

MappedSimpleRead::MappedSimpleRead(MappedSimpleRead&&) = default;

MappedSimpleRead& MappedSimpleRead::operator=(const MappedSimpleRead&) = default;

MappedSimpleRead& MappedSimpleRead::operator=(MappedSimpleRead&&) = default;

MappedSimpleRead::~MappedSimpleRead() = default;

//
// Clipping helpers
//

template <typename T>
T clipContainer(const T& input, const size_t pos, const size_t len)
{
    if (input.empty()) return {};
    return T{input.cbegin() + pos, input.cbegin() + pos + len};
}

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
    const auto clipFrom = result.clipOffset_;
    const auto clipLength = (end - start);
    read.Sequence = clipContainer(read.Sequence, clipFrom, clipLength);
    read.Qualities = clipContainer(read.Qualities, clipFrom, clipLength);
    read.QueryStart = result.qStart_;
    read.QueryEnd = result.qEnd_;
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
    const auto clipFrom = result.clipOffset_;
    const auto clipLength = (end - start);
    read.Sequence = clipContainer(read.Sequence, clipFrom, clipLength);
    read.Qualities = clipContainer(read.Qualities, clipFrom, clipLength);
    read.QueryStart = result.qStart_;
    read.QueryEnd = result.qEnd_;
    read.Cigar = std::move(result.cigar_);
    read.TemplateStart = result.refPos_;
    read.TemplateEnd = read.TemplateStart + ReferenceLength(read.Cigar);
}

void ClipToReference(MappedSimpleRead& read, Position start, Position end,
                     bool exciseFlankingInserts)
{
    // skip out if clip not needed
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
    const auto clipFrom = result.clipOffset_;
    const auto clipLength = (end - start);
    read.Sequence = clipContainer(read.Sequence, clipFrom, clipLength);
    read.Qualities = clipContainer(read.Qualities, clipFrom, clipLength);
    read.QueryStart = result.qStart_;
    read.QueryEnd = result.qEnd_;
    read.Cigar = std::move(result.cigar_);
    read.TemplateStart = result.refPos_;
    read.TemplateEnd = read.TemplateStart + ReferenceLength(read.Cigar);
}

}  // namespace BAM
}  // namespace PacBio
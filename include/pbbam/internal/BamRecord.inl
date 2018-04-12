// File Description
/// \file BamRecord.inl
/// \brief Inline implementations for the BamRecord class.
//
// Author: Derek Barnett

#include "pbbam/BamRecord.h"

namespace PacBio {
namespace BAM {

inline BamRecord BamRecord::Clipped(const BamRecord& input,
                                    const ClipType clipType,
                                    const PacBio::BAM::Position start,
                                    const PacBio::BAM::Position end)
{
    return input.Clipped(clipType, start, end);
}

inline BamRecord BamRecord::Clipped(const ClipType clipType,
                                    const PacBio::BAM::Position start,
                                    const PacBio::BAM::Position end) const
{
    BamRecord result(*this);
    result.Clip(clipType, start, end);
    return result;
}

inline BamRecord BamRecord::Mapped(const BamRecord& input,
                                   const int32_t referenceId,
                                   const Position refStart,
                                   const Strand strand,
                                   const Cigar& cigar,
                                   const uint8_t mappingQuality)
{
    return input.Mapped(referenceId, refStart, strand, cigar, mappingQuality);
}

inline BamRecord BamRecord::Mapped(const int32_t referenceId,
                                   const Position refStart,
                                   const Strand strand,
                                   const Cigar& cigar,
                                   const uint8_t mappingQuality) const
{
    BamRecord result(*this);
    result.Map(referenceId, refStart, strand, cigar, mappingQuality);
    return result;
}

} // namespace BAM
} // namespace PacBio

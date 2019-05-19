// Author: Derek Barnett

#ifndef CLIPPING_H
#define CLIPPING_H

#include "pbbam/Config.h"

#include <pbbam/Cigar.h>
#include <pbbam/Position.h>
#include <pbbam/Strand.h>

namespace PacBio {
namespace BAM {
namespace internal {

struct ClipToQueryConfig
{
    // all clipping
    size_t seqLength_;
    Position original_qStart_;
    Position original_qEnd_;
    Position target_qStart_;
    Position target_qEnd_;

    // for clipping mapped reads
    Position original_tStart_;
    Strand strand_;
    Cigar cigar_;
    bool isMapped_;
};

struct ClipToReferenceConfig : public ClipToQueryConfig
{
    ClipToReferenceConfig(const ClipToQueryConfig& queryConfig, Position originalTEnd,
                          Position targetTStart, Position targetTEnd, bool exciseFlankingInserts);

    Position original_tEnd_;
    Position target_tStart_;
    Position target_tEnd_;
    bool exciseFlankingInserts_;
};

struct ClipResult
{
    ClipResult(size_t clipOffset, Position qStart, Position qEnd);
    ClipResult(size_t clipOffset, Position qStart, Position qEnd, Position refPos, Cigar cigar);

    ClipResult(const ClipResult&);
    ClipResult(ClipResult&&) noexcept;
    ClipResult& operator=(const ClipResult&);
    ClipResult& operator=(ClipResult&&) noexcept;
    ~ClipResult();

    size_t clipOffset_;
    Position qStart_;
    Position qEnd_;
    Position refPos_;
    Cigar cigar_;
};

// configs are non-const so we can steal the input CIGAR, rather than copy,
// but they're otherwise const
ClipResult ClipToQuery(ClipToQueryConfig& config);
ClipResult ClipToReference(ClipToReferenceConfig& config);

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // CLIPPING_H

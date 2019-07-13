// Author: Derek Barnett

#ifndef CLIPPING_H
#define CLIPPING_H

#include "pbbam/Config.h"

#include <pbcopper/data/Cigar.h>
#include <pbcopper/data/Position.h>
#include <pbcopper/data/Strand.h>

namespace PacBio {
namespace BAM {
namespace internal {

struct ClipToQueryConfig
{
    // all clipping
    size_t seqLength_;
    Data::Position original_qStart_;
    Data::Position original_qEnd_;
    Data::Position target_qStart_;
    Data::Position target_qEnd_;

    // for clipping mapped reads
    Data::Position original_tStart_;
    Data::Strand strand_;
    Data::Cigar cigar_;
    bool isMapped_;
};

struct ClipToReferenceConfig : public ClipToQueryConfig
{
    ClipToReferenceConfig(const ClipToQueryConfig& queryConfig, Data::Position originalTEnd,
                          Data::Position targetTStart, Data::Position targetTEnd,
                          bool exciseFlankingInserts);

    Data::Position original_tEnd_;
    Data::Position target_tStart_;
    Data::Position target_tEnd_;
    bool exciseFlankingInserts_;
};

struct ClipResult
{
    ClipResult(size_t clipOffset, Data::Position qStart, Data::Position qEnd);
    ClipResult(size_t clipOffset, Data::Position qStart, Data::Position qEnd, Data::Position refPos,
               Data::Cigar cigar);

    ClipResult(const ClipResult&);
    ClipResult(ClipResult&&) noexcept;
    ClipResult& operator=(const ClipResult&);
    ClipResult& operator=(ClipResult&&) noexcept;
    ~ClipResult();

    size_t clipOffset_;
    Data::Position qStart_;
    Data::Position qEnd_;
    Data::Position refPos_;
    Data::Cigar cigar_;
};

// configs are non-const so we can steal the input CIGAR, rather than copy,
// but they're otherwise const
ClipResult ClipToQuery(ClipToQueryConfig& config);
ClipResult ClipToReference(ClipToReferenceConfig& config);

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // CLIPPING_H


#include "PbbamInternalConfig.h"

#include <cassert>

#include <algorithm>
#include <iostream>

#include <htslib/sam.h>

#include "Clipping.h"

namespace PacBio {
namespace BAM {
namespace internal {

namespace {

// returns reference positions removed from beginning
size_t ClipToQueryImpl(Data::Cigar& cigar, size_t startOffset, size_t endOffset)
{
    size_t refPosRemoved = 0;
    size_t remaining = startOffset;

    // clip CIGAR ops from beginning of query sequence
    while ((remaining > 0) && !cigar.empty()) {
        Data::CigarOperation& op = cigar.front();
        const auto opLength = op.Length();
        const bool consumesQuery = ConsumesQuery(op.Type());
        const bool consumesRef = ConsumesReference(op.Type());

        if (opLength <= remaining) {
            cigar.erase(cigar.begin());
            if (consumesQuery) remaining -= opLength;
            if (consumesRef) refPosRemoved += opLength;
        } else {
            op.Length(opLength - remaining);
            if (consumesRef) refPosRemoved += remaining;
            remaining = 0;
        }
    }

    // clip CIGAR ops from end of query sequence
    remaining = endOffset;
    while ((remaining > 0) && !cigar.empty()) {
        Data::CigarOperation& op = cigar.back();
        const auto opLength = op.Length();
        const bool consumesQuery = ConsumesQuery(op.Type());

        if (opLength <= remaining) {
            cigar.pop_back();
            if (consumesQuery) remaining -= opLength;
        } else {
            op.Length(opLength - remaining);
            remaining = 0;
        }
    }

    return refPosRemoved;
}

void ClipToReferenceImpl(ClipToReferenceConfig& config, size_t* queryPosRemovedFront,
                         size_t* queryPosRemovedBack)
{

    const Data::Position newTStart = std::max(config.original_tStart_, config.target_tStart_);
    const Data::Position newTEnd = std::min(config.original_tEnd_, config.target_tEnd_);

    // fetch a 'working copy' of CIGAR data
    Data::Cigar& cigar = config.cigar_;

    // ------------------------
    // clip leading CIGAR ops
    // ------------------------

    size_t remaining = newTStart - config.original_tStart_;
    while (remaining > 0 && !cigar.empty()) {
        Data::CigarOperation& firstOp = cigar.front();
        const auto firstOpLength = firstOp.Length();
        const bool consumesQuery = ConsumesQuery(firstOp.Type());
        const bool consumesRef = ConsumesReference(firstOp.Type());

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.erase(cigar.begin());
            if (consumesQuery) *queryPosRemovedFront += firstOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or before clip
            if (firstOpLength <= remaining) {
                cigar.erase(cigar.begin());
                if (consumesQuery) *queryPosRemovedFront += firstOpLength;
                if (consumesRef) remaining -= firstOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(firstOpLength > remaining);
                firstOp.Length(firstOpLength - remaining);
                if (consumesQuery) *queryPosRemovedFront += remaining;
                remaining = 0;
            }
        }
    }

    // -------------------------
    // clip trailing CIGAR ops
    // -------------------------

    remaining = config.original_tEnd_ - newTEnd;
    while (remaining > 0 && !cigar.empty()) {
        Data::CigarOperation& lastOp = cigar.back();
        const auto lastOpLength = lastOp.Length();
        const bool consumesQuery = ConsumesQuery(lastOp.Type());
        const bool consumesRef = ConsumesReference(lastOp.Type());

        if (!consumesRef) {

            // e.g. softclip - just pop it completely
            cigar.pop_back();
            if (consumesQuery) *queryPosRemovedBack += lastOpLength;

        } else {
            assert(consumesRef);

            // CIGAR ends at or after clip
            if (lastOpLength <= remaining) {
                cigar.pop_back();
                if (consumesQuery) *queryPosRemovedBack += lastOpLength;
                if (consumesRef) remaining -= lastOpLength;
            }

            // CIGAR straddles clip
            else {
                assert(lastOpLength > remaining);
                lastOp.Length(lastOpLength - remaining);
                if (consumesQuery) *queryPosRemovedBack += remaining;
                remaining = 0;
            }
        }
    }

    if (config.exciseFlankingInserts_) {
        // check for leading insertion
        if (!cigar.empty()) {
            const Data::CigarOperation& op = cigar.front();
            if (op.Type() == Data::CigarOperationType::INSERTION) {
                *queryPosRemovedFront += op.Length();
                cigar.erase(cigar.begin());
            }
        }

        // check for trailing insertion
        if (!cigar.empty()) {
            const Data::CigarOperation& op = cigar.back();
            if (op.Type() == Data::CigarOperationType::INSERTION) {
                *queryPosRemovedBack += op.Length();
                cigar.pop_back();
            }
        }
    }
}

}  // namespace

ClipResult::ClipResult(size_t clipOffset, Data::Position qStart, Data::Position qEnd)
    : clipOffset_{clipOffset}, qStart_{qStart}, qEnd_{qEnd}
{
}

ClipResult::ClipResult(size_t clipOffset, Data::Position qStart, Data::Position qEnd,
                       Data::Position refPos, Data::Cigar cigar)
    : clipOffset_{clipOffset}
    , qStart_{qStart}
    , qEnd_{qEnd}
    , refPos_{refPos}
    , cigar_{std::move(cigar)}
{
}

ClipResult::ClipResult(const ClipResult&) = default;

ClipResult::ClipResult(ClipResult&&) noexcept = default;

ClipResult& ClipResult::operator=(const ClipResult&) = default;

ClipResult& ClipResult::operator=(ClipResult&&) noexcept = default;

ClipResult::~ClipResult() = default;

ClipToReferenceConfig::ClipToReferenceConfig(const ClipToQueryConfig& queryConfig,
                                             Data::Position originalTEnd,
                                             Data::Position targetTStart, Data::Position targetTEnd,
                                             bool exciseFlankingInserts)
    : ClipToQueryConfig{queryConfig}
    , original_tEnd_{originalTEnd}
    , target_tStart_{targetTStart}
    , target_tEnd_{targetTEnd}
    , exciseFlankingInserts_{exciseFlankingInserts}
{
}

ClipResult ClipToQuery(ClipToQueryConfig& config)
{
    // easy out for unmapped reads
    const size_t startOffset = (config.target_qStart_ - config.original_qStart_);
    if (!config.isMapped_)
        return ClipResult{startOffset, config.target_qStart_, config.target_qEnd_};

    // fetch CIGAR (in query orientation)
    Data::Cigar cigar = std::move(config.cigar_);
    if (config.strand_ == Data::Strand::REVERSE) {
        std::reverse(cigar.begin(), cigar.end());
    }

    // do main clipping
    const size_t endOffset = config.original_qEnd_ - config.target_qEnd_;
    const size_t refPosRemoved = ClipToQueryImpl(cigar, startOffset, endOffset);

    // maybe restore CIGAR
    if (config.strand_ == Data::Strand::REVERSE) {
        std::reverse(cigar.begin(), cigar.end());
    }

    // return result
    const Data::Position newPosition = (config.original_tStart_ + refPosRemoved);
    return ClipResult{startOffset, config.target_qStart_, config.target_qEnd_, newPosition,
                      std::move(cigar)};
}

ClipResult ClipToReference(ClipToReferenceConfig& config)
{
    assert(config.isMapped_);

    size_t queryPosRemovedFront = 0;
    size_t queryPosRemovedBack = 0;
    if (config.strand_ == Data::Strand::FORWARD)
        ClipToReferenceImpl(config, &queryPosRemovedFront, &queryPosRemovedBack);
    else
        ClipToReferenceImpl(config, &queryPosRemovedBack, &queryPosRemovedFront);

    const size_t clipOffset = queryPosRemovedFront;
    const Data::Position qStart = config.original_qStart_ + queryPosRemovedFront;
    const Data::Position qEnd = config.original_qEnd_ - queryPosRemovedBack;
    const Data::Position newPos = std::max(config.original_tStart_, config.target_tStart_);
    return ClipResult{clipOffset, qStart, qEnd, newPos, config.cigar_};
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

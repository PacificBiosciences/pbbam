// File Description
/// \file BamRecordView.inl
/// \brief Inline implementations for the BamRecordView class.
//
// Author: Derek Barnett

#include "pbbam/BamRecordView.h"

namespace PacBio {
namespace BAM {

inline BamRecordView::BamRecordView(const BamRecord& record,
                                    const Orientation orientation,
                                    const bool aligned,
                                    const bool exciseSoftClips,
                                    const PulseBehavior pulseBehavior)
    : record_(record)
    , orientation_{orientation}
    , aligned_{aligned}
    , exciseSoftClips_{exciseSoftClips}
    , pulseBehavior_{pulseBehavior}
{ }

inline QualityValues BamRecordView::AltLabelQVs() const
{ return record_.AltLabelQV(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::string BamRecordView::AltLabelTags() const
{ return record_.AltLabelTag(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline QualityValues BamRecordView::DeletionQVs() const
{ return record_.DeletionQV(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::DeletionTags() const
{ return record_.DeletionTag(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::InsertionQVs() const
{ return record_.InsertionQV(orientation_, aligned_, exciseSoftClips_); }

inline Frames BamRecordView::IPD() const
{ return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

inline Frames BamRecordView::PrebaseFrames() const
{ return record_.IPD(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::LabelQVs() const
{ return record_.LabelQV(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline QualityValues BamRecordView::MergeQVs() const
{ return record_.MergeQV(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::PulseMergeQVs() const
{ return record_.PulseMergeQV(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmean() const
{ return record_.Pkmean(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmid() const
{ return record_.Pkmid(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmean2() const
{ return record_.Pkmean2(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::vector<float> BamRecordView::Pkmid2() const
{ return record_.Pkmid2(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline Frames BamRecordView::PrePulseFrames() const
{ return record_.PrePulseFrames(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline std::string BamRecordView::PulseCalls() const
{ return record_.PulseCall(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline Frames BamRecordView::PulseCallWidth() const
{ return record_.PulseCallWidth(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline Frames BamRecordView::PulseWidths() const
{ return record_.PulseWidth(orientation_, aligned_, exciseSoftClips_); }

inline QualityValues BamRecordView::Qualities() const
{ return record_.Qualities(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::Sequence() const
{ return record_.Sequence(orientation_, aligned_, exciseSoftClips_); }

inline std::vector<uint32_t> BamRecordView::StartFrames() const
{ return record_.StartFrame(orientation_, aligned_, exciseSoftClips_, pulseBehavior_); }

inline QualityValues BamRecordView::SubstitutionQVs() const
{ return record_.SubstitutionQV(orientation_, aligned_, exciseSoftClips_); }

inline std::string BamRecordView::SubstitutionTags() const
{ return record_.SubstitutionTag(orientation_, aligned_, exciseSoftClips_); }

} // namespace BAM
} // namespace PacBio

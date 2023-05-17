#ifndef PBBAM_BAMRECORDVIEW_H
#define PBBAM_BAMRECORDVIEW_H

#include <pbbam/Config.h>

#include <pbbam/BamRecord.h>

#include <cstdint>

namespace PacBio {
namespace BAM {

/// \brief Provides a re-usable "view" onto a BamRecord
///
/// This class acts a convenience wrapper for working with per-base BamRecord
/// data. Most of these BamRecord methods take a list of parameters, to adjust
/// how the underlying data are presented to client code. Often these parameters
/// will be re-used for each BamRecord method call. Thus, to simplify such
/// client code, a BamRecordView can be used to state those parameters once, and
/// then simply request the desired fields.
///
/// \internal
/// \todo Sync up method names with BamRecord
/// \endinternal
///
class PBBAM_EXPORT BamRecordView
{
public:
    /// \brief Constructs a view onto \p record using the supplied parameters.
    ///
    /// For frame or QV data, if \p aligned is true, a value of 0 (Accuracy or
    /// QualityValue) will be used at each inserted or padded base location.
    ///
    /// \param[in] record           BamRecord data source.
    /// \param[in] orientation      Orientation of output.
    /// \param[in] aligned          if true, gaps/padding will be inserted, per
    ///                             Cigar info.
    /// \param[in] exciseSoftClips  if true, any soft-clipped positions will be
    ///                             removed from query ends
    ///
    BamRecordView(const BamRecord& record, Data::Orientation orientation, bool aligned,
                  bool exciseSoftClips, PulseBehavior pulseBehavior = PulseBehavior::ALL);

public:
    /// \returns BamRecord::AltLabelQV with this view's parameters applied
    Data::QualityValues AltLabelQVs() const;

    /// \returns BamRecord::AltLabelTag with this view's parameters applied
    std::string AltLabelTags() const;

    /// \returns BamRecord::DeletionQV with this view's parameters applied
    Data::QualityValues DeletionQVs() const;

    /// \returns BamRecord::DeletionTag with this view's parameters applied
    std::string DeletionTags() const;

    /// \returns BamRecord::InsertionQV with this view's parameters applied
    Data::QualityValues InsertionQVs() const;

    /// \returns BamRecord::IPD with this view's parameters applied
    Data::Frames IPD() const;

    /// \returns BamRecord::LabelQV with this view's parameters applied
    Data::QualityValues LabelQVs() const;

    /// \returns BamRecord::MergeQV with this view's parameters applied
    Data::QualityValues MergeQVs() const;

    /// \returns BamRecord::PulseMergeQV with this view's parameters applied
    Data::QualityValues PulseMergeQVs() const;

    /// \returns BamRecord::Pkmean with this view's parameters applied
    std::vector<float> Pkmean() const;

    /// \returns BamRecord::Pkmid with this view's parameters applied
    std::vector<float> Pkmid() const;

    /// \returns BamRecord::Pkmean2 with this view's parameters applied
    std::vector<float> Pkmean2() const;

    /// \returns BamRecord::Pkmid2 with this view's parameters applied
    std::vector<float> Pkmid2() const;

    /// \returns BamRecord::PreBaseFrames with this view's parameters applied
    Data::Frames PrebaseFrames() const;

    /// \returns BamRecord::PrePulseFrames with this view's parameters applied
    Data::Frames PrePulseFrames() const;

    /// \returns BamRecord::PulseCalls with this view's parameters applied
    std::string PulseCalls() const;

    /// \returns BamRecord::PulseCallWidth with this view's parameters applied
    Data::Frames PulseCallWidth() const;

    /// \returns BamRecord::PulseWidths with this view's parameters applied
    Data::Frames PulseWidths() const;

    /// \returns BamRecord::Qualities with this view's parameters applied
    Data::QualityValues Qualities() const;

    /// \returns BamRecord::Sequence with this view's parameters applied
    std::string Sequence() const;

    /// \returns BamRecord::StartFrame with this view's parameters applied
    std::vector<std::uint32_t> StartFrames() const;

    /// \returns BamRecord::SubstitutionQV with this view's parameters applied
    Data::QualityValues SubstitutionQVs() const;

    /// \returns BamRecord::SubstitutionTag with this view's parameters applied
    std::string SubstitutionTags() const;

private:
    const BamRecord& record_;
    Data::Orientation orientation_;
    bool aligned_;
    bool exciseSoftClips_;
    PulseBehavior pulseBehavior_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_BAMRECORDVIEW_H

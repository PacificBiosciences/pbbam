// File Description
/// \file FastqSequence.inl
/// \brief Inline implementations for the FastqSequence class.
//
// Author: Derek Barnett

#include "pbbam/FastqSequence.h"

namespace PacBio {
namespace BAM {

inline FastqSequence::FastqSequence(std::string name,
                                    std::string bases,
                                    QualityValues qualities)
    : FastaSequence{std::move(name), std::move(bases)}
    , qualities_{std::move(qualities)}
{ }

inline FastqSequence::FastqSequence(std::string name,
                                    std::string bases,
                                    std::string qualities)
    : FastaSequence{std::move(name), std::move(bases)}
    , qualities_{QualityValues::FromFastq(qualities)}
{ }

inline const QualityValues& FastqSequence::Qualities() const
{ return qualities_; }

} // namespace BAM
} // namespace PacBio

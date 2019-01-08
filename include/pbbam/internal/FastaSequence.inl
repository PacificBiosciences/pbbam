// File Description
/// \file FastaSequence.inl
/// \brief Inline implementations for the FastaSequence class.
//
// Author: Derek Barnett

#include "pbbam/FastaSequence.h"

namespace PacBio {
namespace BAM {

inline FastaSequence::FastaSequence(std::string name,
                                    std::string bases)
    : name_{std::move(name)}
    , bases_{std::move(bases)}
{ }

inline const std::string& FastaSequence::Bases() const
{ return bases_; }

inline const std::string& FastaSequence::Name() const
{ return name_; }

} // namespace BAM
} // namespace PacBio

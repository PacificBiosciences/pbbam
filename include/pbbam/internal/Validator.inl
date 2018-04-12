// File Description
/// \file Validator.inl
/// \brief Inline implementations for the Validator class.
//
// Author: Derek Barnett

#include "pbbam/Validator.h"
#include <stdexcept>

namespace PacBio {
namespace BAM {

inline bool Validator::IsValid(const BamFile& file, const bool entireFile)
{
    try {
        if (entireFile)
            ValidateEntireFile(file, 1);
        else
            ValidateFileMetadata(file, 1);
        return true;
    } catch (std::exception&) {
        return false;
    }
}

inline bool Validator::IsValid(const BamHeader& header)
{
    try {
        Validate(header, 1);
        return true;
    } catch (std::exception&) {
        return false;
    }
}

inline bool Validator::IsValid(const BamRecord& record)
{
    try {
        Validate(record, 1);
        return true;
    } catch (std::exception&) {
        return false;
    }
}

inline bool Validator::IsValid(const ReadGroupInfo& rg)
{
    try {
        Validate(rg, 1);
        return true;
    } catch (std::exception&) {
        return false;
    }
}

} // namespace BAM
} // namespace PacBio

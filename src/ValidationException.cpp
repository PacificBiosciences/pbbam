// File Description
/// \file ValidationException.cpp
/// \brief Implements the ValidationException class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/exception/ValidationException.h"

namespace PacBio {
namespace BAM {

ValidationException::ValidationException(ErrorMap fileErrors, ErrorMap readGroupErrors,
                                         ErrorMap recordErrors)
    : std::runtime_error{""}
    , fileErrors_{std::move(fileErrors)}
    , readGroupErrors_{std::move(readGroupErrors)}
    , recordErrors_{std::move(recordErrors)}
{
    FormatMessage();
}

const ValidationException::ErrorMap& ValidationException::FileErrors() const { return fileErrors_; }

const ValidationException::ErrorMap& ValidationException::ReadGroupErrors() const
{
    return readGroupErrors_;
}

const ValidationException::ErrorMap& ValidationException::RecordErrors() const
{
    return recordErrors_;
}

const char* ValidationException::what() const noexcept { return msg_.c_str(); }

void ValidationException::FormatMessage()
{
    std::ostringstream s;
    s << "Validation failed:\n";

    // file errors
    if (!fileErrors_.empty()) {
        for (const auto& fileError : fileErrors_) {
            s << "  In file (" << fileError.first << ") : \n";
            for (const auto& e : fileError.second)
                s << "    " << e << '\n';
        }
    }

    // read group errors
    if (!readGroupErrors_.empty()) {
        for (const auto& rgError : readGroupErrors_) {
            s << "  In read group (" << rgError.first << ") :\n";
            for (const auto& e : rgError.second)
                s << "    " << e << '\n';
        }
    }

    // record errors
    if (!recordErrors_.empty()) {
        for (const auto& recordError : readGroupErrors_) {
            s << "  In record (" << recordError.first << ") : \n";
            for (const auto& e : recordError.second)
                s << "    " << e << '\n';
        }
    }

    msg_ = s.str();
}

}  // namespace BAM
}  // namespace PacBio

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
    : std::runtime_error("")
    , fileErrors_(std::move(fileErrors))
    , readGroupErrors_(std::move(readGroupErrors))
    , recordErrors_(std::move(recordErrors))
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
    std::stringstream s;
    s << "Validation failed: " << std::endl;

    // file errors
    if (!fileErrors_.empty()) {
        auto fileIter = fileErrors_.cbegin();
        auto fileEnd = fileErrors_.cend();
        for (; fileIter != fileEnd; ++fileIter) {
            s << "  In file (" << fileIter->first << ") : " << std::endl;
            const auto& errors = fileIter->second;
            for (const auto& e : errors)
                s << "    " << e << std::endl;
        }
    }

    // read group errors
    if (!readGroupErrors_.empty()) {
        auto rgIter = readGroupErrors_.cbegin();
        auto rgEnd = readGroupErrors_.cend();
        for (; rgIter != rgEnd; ++rgIter) {
            s << "  In read group (" << rgIter->first << ") : " << std::endl;
            const auto& errors = rgIter->second;
            for (const auto& e : errors)
                s << "    " << e << std::endl;
        }
    }

    // record errors
    if (!recordErrors_.empty()) {
        auto recIter = recordErrors_.cbegin();
        auto recEnd = recordErrors_.cend();
        for (; recIter != recEnd; ++recIter) {
            s << "  In record (" << recIter->first << ") : " << std::endl;
            const auto& errors = recIter->second;
            for (const auto& e : errors)
                s << "    " << e << std::endl;
        }
    }

    msg_ = s.str();
}

}  // namespace BAM
}  // namespace PacBio

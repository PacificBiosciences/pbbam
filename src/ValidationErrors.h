// File Description
/// \file ValidationErrors.h
/// \brief Defines the ValidationErrors class.
//
// Author: Derek Barnett

#ifndef VALIDATIONERRORS_H
#define VALIDATIONERRORS_H

#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {
namespace internal {

/// The ValidationErrors class catches error messages accumulated during
/// validation (see Validator).
///
/// Convenience methods are provided for different BAM components, to help
/// format the displayed output.
///
/// A maximum number of errors can be provided at construction, and this class
/// will automatially throw a ValidationException whenever that count is reached.
/// Otherwise, the Validator will check IsEmpty() and call ThrowErrors() if true.
///
class ValidationErrors
{
public:
    typedef std::vector<std::string> ErrorList;
    typedef std::map<std::string, ErrorList> ErrorMap;

public:
    static const size_t MAX = std::numeric_limits<size_t>::max();

public:
    ValidationErrors(const size_t maxNumErrors = ValidationErrors::MAX);

public:
    void AddFileError(const std::string& fn, std::string details);
    void AddReadGroupError(const std::string& rg, std::string details);
    void AddRecordError(const std::string& name, std::string details);
    void AddTagLengthError(const std::string& name, const std::string& tagLabel,
                           const std::string& tagName, const size_t observed,
                           const size_t expected);

public:
    bool IsEmpty() const;
    size_t MaxNumErrors() const;
    void ThrowErrors();

private:
    size_t maxNumErrors_;
    size_t currentNumErrors_;
    ErrorMap fileErrors_;
    ErrorMap readGroupErrors_;
    ErrorMap recordErrors_;

private:
    void OnErrorAdded();
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // VALIDATIONERRORS_H

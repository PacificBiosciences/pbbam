// File Description
/// \file ValidationException.h
/// \brief Defines the ValidationException class.
//
// Author: Derek Barnett

#ifndef VALIDATIONEXCEPTION_H
#define VALIDATIONEXCEPTION_H

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief The ValidationExecption represents an exception that will be thrown
///        when any error is encountered using the Validator API. In addition to
///        a default display message, it provides programmatic access to all
///        reported error messages.
///
/// \sa Validator::Validate(const BamRecord& record)
///
class ValidationException : public std::runtime_error
{
public:
    using ErrorList = std::vector<std::string>;
    using ErrorMap = std::map<std::string, ErrorList>;

public:
    ValidationException(ErrorMap fileErrors, ErrorMap readGroupErrors, ErrorMap recordErrors);

    // This is a work around for the Intel PHI compiler (icpc)
    ~ValidationException() throw() {}

public:
    const ErrorMap& FileErrors() const;
    const ErrorMap& ReadGroupErrors() const;
    const ErrorMap& RecordErrors() const;

    const char* what() const noexcept override;

private:
    ErrorMap fileErrors_;
    ErrorMap readGroupErrors_;
    ErrorMap recordErrors_;
    std::string msg_;

private:
    void FormatMessage();
};

}  // namespace BAM
}  // namespace PacBio

#endif  // VALIDATIONEXCEPTION_H

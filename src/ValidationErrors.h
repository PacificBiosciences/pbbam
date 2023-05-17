#ifndef PBBAM_VALIDATIONERRORS_H
#define PBBAM_VALIDATIONERRORS_H

#include <pbbam/Config.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace BAM {

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

    static const std::size_t MAX = std::numeric_limits<std::size_t>::max();

    explicit ValidationErrors(std::size_t maxNumErrors = ValidationErrors::MAX);

    void AddFileError(const std::string& fn, std::string details);
    void AddReadGroupError(const std::string& rg, std::string details);
    void AddRecordError(const std::string& name, std::string details);
    void AddTagLengthError(const std::string& name, const std::string& tagLabel,
                           const std::string& tagName, std::size_t observed, std::size_t expected);

    bool IsEmpty() const;
    std::size_t MaxNumErrors() const;
    void ThrowErrors();

private:
    std::size_t maxNumErrors_;
    std::size_t currentNumErrors_;
    ErrorMap fileErrors_;
    ErrorMap readGroupErrors_;
    ErrorMap recordErrors_;

    void OnErrorAdded();
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_VALIDATIONERRORS_H

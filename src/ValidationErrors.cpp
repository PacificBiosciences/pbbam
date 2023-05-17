#include "PbbamInternalConfig.h"

#include <pbbam/exception/ValidationException.h>

#include <pbbam/StringUtilities.h>
#include "ValidationErrors.h"

#include <sstream>

#include <cstddef>

namespace PacBio {
namespace BAM {

const std::size_t ValidationErrors::MAX;

ValidationErrors::ValidationErrors(const std::size_t maxNumErrors)
    : maxNumErrors_{maxNumErrors}, currentNumErrors_{0}
{
    if (maxNumErrors_ == 0) {
        maxNumErrors_ = ValidationErrors::MAX;
    }
}

void ValidationErrors::AddFileError(const std::string& fn, std::string details)
{
    fileErrors_[fn].push_back(std::move(details));
    OnErrorAdded();
}

void ValidationErrors::AddReadGroupError(const std::string& rg, std::string details)
{
    readGroupErrors_[rg].push_back(std::move(details));
    OnErrorAdded();
}

void ValidationErrors::AddRecordError(const std::string& name, std::string details)
{
    recordErrors_[name].push_back(std::move(details));
    OnErrorAdded();
}

void ValidationErrors::AddTagLengthError(const std::string& name, const std::string& tagLabel,
                                         const std::string& tagName, const std::size_t observed,
                                         const std::size_t expected)
{
    // format
    std::ostringstream s;
    s << tagLabel << " tag (" << tagName << ") length: " << observed
      << ", does not match expected length: " << expected;
    AddRecordError(name, s.str());
}

bool ValidationErrors::IsEmpty() const { return currentNumErrors_ == 0; }

size_t ValidationErrors::MaxNumErrors() const { return maxNumErrors_; }

void ValidationErrors::OnErrorAdded()
{
    ++currentNumErrors_;
    if (currentNumErrors_ == maxNumErrors_) {
        ThrowErrors();
    }
}

void ValidationErrors::ThrowErrors()
{
    throw ValidationException{std::move(fileErrors_), std::move(readGroupErrors_),
                              std::move(recordErrors_)};
}

}  // namespace BAM
}  // namespace PacBio

// File Description
/// \file InvalidSequencingChemistryException.h
/// \brief Defines the InvalidSequencingChemistryException class.
//
// Author: Derek Barnett

#ifndef INVALIDSEQUENCINGCHEMISTRYEXCEPTION_H
#define INVALIDSEQUENCINGCHEMISTRYEXCEPTION_H

#include <exception>
#include <sstream>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The InvalidSequencingChemistryException class represents an exception
///        that will be thrown when an invalid sequencing chemistry combination
///        is encountered.
///
class InvalidSequencingChemistryException : public std::exception
{
public:
    InvalidSequencingChemistryException(std::string bindingKit, std::string sequencingKit,
                                        std::string basecallerVersion)
        : bindingKit_(std::move(bindingKit))
        , sequencingKit_(std::move(sequencingKit))
        , basecallerVersion_(std::move(basecallerVersion))
    {
        std::ostringstream s;
        s << "unsupported sequencing chemistry combination:\n"
          << "    binding kit:        " << bindingKit_ << '\n'
          << "    sequencing kit:     " << sequencingKit_ << '\n'
          << "    basecaller version: " << basecallerVersion_ << '\n';
        what_ = s.str();
    }

    // This is a work around for the Intel PHI compiler (icpc)
    ~InvalidSequencingChemistryException() throw() {}

public:
    const std::string& BindingKit() const { return bindingKit_; }

    const std::string& SequencingKit() const { return sequencingKit_; }

    const std::string& BasecallerVersion() const { return basecallerVersion_; }

public:
    const char* what() const noexcept override { return what_.c_str(); }

protected:
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    std::string what_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // INVALIDSEQUENCINGCHEMISTRYEXCEPTION_H

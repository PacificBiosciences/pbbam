// File Description
/// \file BundleChemistryMappingException.h
/// \brief Defines the BundleChemistryMappingException class.
//
// Author: Derek Barnett, Lance Hepler

#ifndef BUNDLECHEMISTRYMAPPINGEXCEPTION_H
#define BUNDLECHEMISTRYMAPPINGEXCEPTION_H

#include <exception>
#include <string>

namespace PacBio {
namespace BAM {

/// \brief The BundleChemistryMappingException class represents an exception
///        that will be thrown when an invalid sequencing chemistry combination
///        is encountered.
///
class BundleChemistryMappingException : public std::exception
{
public:
    BundleChemistryMappingException(std::string mappingXml, std::string msg)
        : mappingXml_(std::move(mappingXml))
        , what_(std::string("invalid ") + mappingXml_ + ": " + std::move(msg))
    {
    }

    // This is a work around for the Intel PHI compiler (icpc)
    ~BundleChemistryMappingException() throw() {}

public:
    const char* what() const noexcept override { return what_.c_str(); }

protected:
    std::string mappingXml_;
    std::string what_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // BUNDLECHEMISTRYMAPPINGEXCEPTION_H

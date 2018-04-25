// File Description
/// \file Validator.h
/// \brief Defines the Validator class.
//
// Author: Derek Barnett

#ifndef VALIDATOR_H
#define VALIDATOR_H

#include <cstddef>
#include <limits>
#include "pbbam/Config.h"
#include "pbbam/exception/ValidationException.h"

namespace PacBio {
namespace BAM {

class BamFile;
class BamHeader;
class BamRecord;
class ReadGroupInfo;

/// \brief The Validator class provides validation for %BAM data.
///
/// There are 2 ways to use this class. If you are only compared with a quick &
/// dirty, yes/no validation, then you can use the IsValid() methods. This will
/// swallow the specific cause of the failure, but you don't have to catch an
/// exception and handle it in your client code. If you want to know,
/// specifically, what failed, then you can use the Validate*() methods that
/// will throw a ValidationException if the object is invalid. This exception
/// will provide more details as to what failed and why.
///
/// See documentation for Config.h for details on building pbbam with
/// auto-validation enabled.
///
class PBBAM_EXPORT Validator
{
public:
    /// \brief Checks that a %BAM file conforms to the %PacBio specification.
    ///
    /// When \p entireFile is false, this method only checks file metadata. If
    /// \p entireFile is true, all records are checked as well.
    ///
    /// \param[in] file         %BAM header to validate
    /// \param[in] entireFile   check records in addition to metadata
    /// \returns true if \p file passes validation checks
    ///
    /// \sa Validator::ValidateFileMetdata, Validator::ValidateEntireFile
    ///
    static bool IsValid(const BamFile& file, const bool entireFile);

    /// \brief Checks that a %BAM header conforms to the %PacBio specification.
    ///
    /// \returns true if \p header passes validation checks
    ///
    /// \sa Validator::Validate(const BamHeader& header)
    ///
    static bool IsValid(const BamHeader& header);

    /// \brief Checks that a %BAM read group conforms to the %PacBio
    ///        specification.
    ///
    /// \returns true if \p rg passes validation checks
    ///
    /// \sa Validator::Validate(const ReadGroupInfo& rg)
    ///
    static bool IsValid(const ReadGroupInfo& rg);

    /// \brief Checks that a %BAM record conforms to the %PacBio specification.
    ///
    /// \returns true if \p record passes validation checks
    ///
    /// \sa Validator::Validate(const BamRecord& record)
    ///
    static bool IsValid(const BamRecord& record);

public:
    Validator() = delete;

    /// \brief Checks that a %BAM file's header conforms to the
    ///        %PacBio specification.
    ///
    /// This validation step checks the SAM/%BAM version number, sort order,
    /// PacBioBAM version number, and calls Validate(readGroup) internally for
    /// all read groups.
    ///
    /// \param[in] file         %BAM header to validate
    /// \param[in] maxErrors    maximum number of errors to allow before throwing
    ///
    /// \throws ValidationException if \p header fails validation checks
    ///
    static void Validate(const BamHeader& header,
                         const size_t maxErrors = std::numeric_limits<size_t>::max());

    /// \brief Checks that a %BAM read group conforms to the %PacBio
    ///        specification.
    ///
    /// \param[in] rg           %BAM read group to validate
    /// \param[in] maxErrors    maximum number of errors to allow before throwing
    ///
    /// \throws ValidationException if \p rg fails validation checks
    ///
    static void Validate(const ReadGroupInfo& rg,
                         const size_t maxErrors = std::numeric_limits<size_t>::max());

    /// \brief Checks that a %BAM record conforms to the %PacBio specification.
    ///
    /// \param[in] record       %BAM record to validate
    /// \param[in] maxErrors    maximum number of errors to allow before throwing
    ///
    /// \throws ValidationException if \p record fails validation checks
    ///
    static void Validate(const BamRecord& record,
                         const size_t maxErrors = std::numeric_limits<size_t>::max());

    /// \brief Checks that a %BAM file's (entire) contents conform to the
    ///        %PacBio specification.
    ///
    /// This is equivalent to:
    ///
    /// \code
    /// Validator::ValidateMetadata(file);
    /// EntireFileQuery query(file);
    /// for (const BamRecord& record : query)
    ///     Validator::Validate(record);
    /// \endcode
    ///
    /// \param[in] file         %BAM file to validate
    /// \param[in] maxErrors    maximum number of errors to allow before throwing
    ///
    /// \throws ValidationException if \p file fails validation checks
    ///
    static void ValidateEntireFile(const BamFile& file,
                                   const size_t maxErrors = std::numeric_limits<size_t>::max());

    /// \brief Checks that a %BAM file's metadata conforms to the
    ///        %PacBio specification.
    ///
    /// This validation step checks the filename, ensures EOF marker, and
    /// presence of PBI. It also calls Validate(file.Header()) internally.
    ///
    /// \param[in] file         %BAM header to validate
    /// \param[in] maxErrors    maximum number of errors to allow before throwing
    ///
    /// \throws ValidationException if \p header fails validation checks
    ///
    static void ValidateFileMetadata(const BamFile& file,
                                     const size_t maxErrors = std::numeric_limits<size_t>::max());
};

}  // namespace BAM
}  // namespace PacBio

#include "internal/Validator.inl"

#endif  // VALIDATOR_H

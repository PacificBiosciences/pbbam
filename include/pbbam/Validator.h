// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.
//
// File Description
/// \file Validator.h
/// \brief Defines the Validator class.
//
// Author: Derek Barnett

#ifndef VALIDATOR_H
#define VALIDATOR_H

#include "pbbam/Config.h"
#include "pbbam/exception/ValidationException.h"
#include <cstddef>
#include <limits>

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

} // namespace BAM
} // namespace PacBio

#include "internal/Validator.inl"

#endif // VALIDATOR_H

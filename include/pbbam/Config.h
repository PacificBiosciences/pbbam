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
/// \file Config.h
/// \brief Defines library-wide macros & global variables.
//
// Author: Derek Barnett

#ifndef PBBAM_CONFIG_H
#define PBBAM_CONFIG_H

#include <cstdint>

#ifndef INT8_MAX
#define INT8_MAX         127
#endif
#ifndef INT16_MAX
#define INT16_MAX        32767
#endif
#ifndef INT32_MAX
#define INT32_MAX        2147483647
#endif
#ifndef INT64_MAX
#define INT64_MAX        9223372036854775807LL
#endif
#ifndef INT8_MIN
#define INT8_MIN          -128
#endif
#ifndef INT16_MIN 
#define INT16_MIN         -32768
#endif
#ifndef INT32_MIN
#define INT32_MIN        (-INT32_MAX-1)
#endif
#ifndef INT64_MIN
#define INT64_MIN        (-INT64_MAX-1)
#endif
#ifndef UINT8_MAX 
#define UINT8_MAX         255
#endif
#ifndef UINT16_MAX
#define UINT16_MAX        65535
#endif
#ifndef UINT32_MAX
#define UINT32_MAX        4294967295U
#endif
#ifndef UINT64_MAX
#define UINT64_MAX        18446744073709551615ULL
#endif

/// \name Library Import/Export
/// \{

#ifndef PBBAM_LIBRARY_EXPORT
#  if defined(WIN32)
#    define PBBAM_LIBRARY_EXPORT __declspec(dllexport)
#  else
#    define PBBAM_LIBRARY_EXPORT __attribute__((visibility("default")))
#  endif
#endif

#ifndef PBBAM_LIBRARY_IMPORT
#  if defined(WIN32)
#    define PBBAM_LIBRARY_IMPORT __declspec(dllimport)
#  else
#    define PBBAM_LIBRARY_IMPORT
#  endif
#endif

#ifndef PBBAM_EXPORT
#  if defined(PBBAM_LIBRARY)
#    define PBBAM_EXPORT PBBAM_LIBRARY_EXPORT
#  else
#    define PBBAM_EXPORT PBBAM_LIBRARY_IMPORT
#  endif
#endif

/// \}

/// \name Shared Pointer Settings
/// \{

// uncomment this define, or pass via command-line (-DPBBAM_USE_BOOST_SHARED_PTR),
// to use boost::shared_ptr<T> instead of std::shared_ptr<T>
//
//#define PBBAM_USE_BOOST_SHARED_PTR

#ifdef PBBAM_USE_BOOST_SHARED_PTR
#  include <boost/smart_ptr/shared_ptr.hpp>
#  define PBBAM_SHARED_PTR boost::shared_ptr
#else
#  include <memory>
#  define PBBAM_SHARED_PTR std::shared_ptr
#endif

/// \}

/// \name Class Definition Helpers
/// \{

/// \brief Disables the use of copy constructors and assignment operators for a
///        class.
///
/// To use, place the macro in a class's private section:
/// \code{.cpp}
/// struct Foo {
/// private:
///     DISABLE_COPY(Foo);
/// };
/// \endcode
///
#ifndef DISABLE_COPY
#define DISABLE_COPY(Class) \
    Class(const Class&); \
    Class& operator=(const Class&)
#endif

/// \brief Disables the use of move constructors and assignment operators for a
///        class.
///
/// To use, place the macro in a class's private section:
/// \code{.cpp}
/// struct Foo {
/// private:
///     DISABLE_MOVE(Foo);
/// };
/// \endcode
///
#ifndef DISABLE_MOVE
#define DISABLE_MOVE(Class) \
    Class(Class&&); \
    Class& operator=(Class&&);
#endif

/// \brief Disables the use of copy & move constructors and assignment operators f
///        or a class.
///
/// To use, place the macro in a class's private section:
/// \code{.cpp}
/// struct Foo {
/// private:
///     DISABLE_MOVE_AND_COPY(Foo);
/// };
/// \endcode
///
#ifndef DISABLE_MOVE_AND_COPY
#define DISABLE_MOVE_AND_COPY(Class) \
    DISABLE_MOVE(Class) \
    DISABLE_COPY(Class)
#endif

/// \}

// \brief Auto-validation
//
// To validate BAM components (header, records, etc.) you can either use the
// Validator API provided, or enable auto-validation. To compile pbbam for
// auto-validation, add the -DPacBioBAM_auto_validate=ON option to your cmake
// invocation.
//
//
#ifndef PBBAM_AUTOVALIDATE
#  define PBBAM_AUTOVALIDATE 0
#endif

/// \}

namespace PacBio {
namespace BAM {

/// \name Verbosity Settings
/// \{

/// \brief Sets the desired verbosity level of htslib warnings.
///
/// Change this value to allow debug/warning statements from htslib itself.
/// The valid range seems to be [0-3], where 0 indicates OFF, and 3 is the
/// most verbose.
///
/// By default, pbbam disables htslib statements to keep output channels clean.
/// We rely on exceptions & their associated messages instead.
///
/// This global variable is obviously not thread-safe by any means. But as a
/// debug flag, it is unlikely to cause any real issues. The worst case would be
/// unexpected presence/absence of output statements.
///
extern int HtslibVerbosity;

/// \}

} // namespace BAM
} // namespace PacBio

#endif // PBBAM_CONFIG_H

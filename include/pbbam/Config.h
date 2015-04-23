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

// Author: Derek Barnett

#ifndef PBBAM_CONFIG_H
#define PBBAM_CONFIG_H

// --------------------------------
// standard types
// --------------------------------

#include <cstdint>

// -------------------------------------
// library symbol import/export macros
// -------------------------------------

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

// ----------------------------------------------------
// setup the shared_ptr implementation we'll be using
// ----------------------------------------------------

// uncomment this define, or pass via command-line (-DPBBAM_USE_BOOST_SHARED_PTR),
// to use boost::shared_ptr<T> instead of std::shared_ptr<T>
//#define PBBAM_USE_BOOST_SHARED_PTR

#ifdef PBBAM_USE_BOOST_SHARED_PTR
#  include <boost/smart_ptr/shared_ptr.hpp>
#  define PBBAM_SHARED_PTR boost::shared_ptr
#else
#  include <memory>
#  define PBBAM_SHARED_PTR std::shared_ptr
#endif

// ----------------------------------------------------
// htslib verbosity level
// ----------------------------------------------------

namespace PacBio {
namespace BAM {

/// \brief Sets the desired verbosity level of htslib warnings.
///
/// Change this value to allow debug/warning statements from htslib.
/// The valid range seems to be [0-3], where 0->OFF, and 3->most verbose.
///
extern int HtslibVerbosity;

} // namespace BAM
} // namespace PacBio

// ----------------------------------------------------
// additional helper macros
// ----------------------------------------------------

#ifndef DISABLE_COPY
#define DISABLE_COPY(Class) \
    Class(const Class&); \
    Class& operator=(const Class&)
#endif

#ifndef DISABLE_MOVE_AND_COPY
#define DISABLE_MOVE_AND_COPY(Class) \
    Class(Class&&); \
    Class& operator=(Class&&); \
    DISABLE_COPY(Class)
#endif

#endif // PBBAM_CONFIG_H

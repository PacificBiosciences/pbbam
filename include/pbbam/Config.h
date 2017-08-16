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

/// \name Library Import/Export
/// \{

#ifndef PBBAM_EXPORT
#  if defined(WIN32)
#    define PBBAM_EXPORT __declspec(dllimport)
#  else
#    define PBBAM_EXPORT
#  endif
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

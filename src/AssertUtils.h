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

#ifndef ASSERTUTILS_H
#define ASSERTUTILS_H

// ---------------------------------------------------
// This file contains dev/debugging helper utilities
// ---------------------------------------------------

#ifndef PBBAM_UNUSED
#  define PBBAM_UNUSED(x) (void)x;
#endif

namespace PacBio {
namespace BAM {
namespace internal {

inline void pbbam_noop(void) { }

// a la fprintf(...). Auto-adds a newline
void printError(const char* msg, ...);
void printInfo(const char* msg, ...);
void printFailedAssert(const char* msg);

} // namespace internal
} // namespace BAM
} // namespace PacBio

//
// This assert construct below allows us to report failures as well as take some
// fallback action (return, break, continue, etc) so as not to crash at runtime.
// In other words, it's basically a 'weak' assert with customized information &
// failure response.
//
// PB_VERIFY(cond)           if condition fails, print message
// PB_ASSERT(cond, action)   if condition fails, print message & perform action
// PB_ASSERT_OR_BREAK        overload of ASSERT where action is 'break'
// PB_ASSERT_OR_CONTINUE     overload of ASSERT where action is 'continue'
// PB_ASSERT_OR_RETURN       overload of ASSERT where action is 'return'
// PB_ASSERT_OR_RETURN_VALUE overload of ASSERT where action is 'return <value>'
// PB_ASSERT_UNREACHABLE     overload of ASSERT(false) where action is a no-op. Used as a visual marker for
//                        unreachable code-paths (e.g. invalid values in a switch statement)
//
#define PB_ASSERT_STRINGIFY2(x) #x
#define PB_ASSERT_STRINGIFY(x) PB_ASSERT_STRINGIFY2(x)
#define PB_ASSERT_STRING(cond) ::PacBio::BAM::internal::printFailedAssert( \
    "\"" cond"\" in file " __FILE__ ", line " PB_ASSERT_STRINGIFY(__LINE__))

#define PB_VERIFY(cond)             if (cond) {} else { PB_ASSERT_STRING(#cond);         } do {} while (0)
#define PB_ASSERT(cond, action)     if (cond) {} else { PB_ASSERT_STRING(#cond); action; } do {} while (0)
#define PB_ASSERT_OR_BREAK(cond)    PB_ASSERT(cond, break)
#define PB_ASSERT_OR_CONTINUE(cond) PB_ASSERT(cond, continue)
#define PB_ASSERT_OR_RETURN(cond)   PB_ASSERT(cond, return)
#define PB_ASSERT_OR_RETURN_VALUE(cond, value) PB_ASSERT(cond, return value)

#define PB_ASSERT_UNREACHABLE PB_ASSERT(false, ::PacBio::BAM::internal::pbbam_noop())

#endif // ASSERTUTILS_H

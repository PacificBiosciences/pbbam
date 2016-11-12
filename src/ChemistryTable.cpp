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

// Author: Lance Hepler

#include "ChemistryTable.h"

namespace PacBio {
namespace BAM {
namespace internal {

extern const std::vector<std::array<std::string, 4>> ChemistryTable = {

    // BindingKit, SequencingKit, BasecallerVersion, Chemistry

    // RS
    {{"100356300",   "100356200",   "2.1", "P6-C4"}},
    {{"100356300",   "100356200",   "2.3", "P6-C4"}},
    {{"100356300",   "100612400",   "2.1", "P6-C4"}},
    {{"100356300",   "100612400",   "2.3", "P6-C4"}},
    {{"100372700",   "100356200",   "2.1", "P6-C4"}},
    {{"100372700",   "100356200",   "2.3", "P6-C4"}},
    {{"100372700",   "100612400",   "2.1", "P6-C4"}},
    {{"100372700",   "100612400",   "2.3", "P6-C4"}},

    // 3.0 ("Dromedary"): S/P1-C1/beta
    {{"100-619-300", "100-620-000", "3.0", "S/P1-C1/beta"}},
    {{"100-619-300", "100-620-000", "3.1", "S/P1-C1/beta"}},

    // 3.1 ("Echidna"): S/P1-C1.1
    {{"100-619-300", "100-867-300", "3.1", "S/P1-C1.1"}},
    {{"100-619-300", "100-867-300", "3.2", "S/P1-C1.1"}},
    {{"100-619-300", "100-867-300", "3.3", "S/P1-C1.1"}},

    // 3.1.1 ("Flea"): S/P1-C1.2
    {{"100-619-300", "100-902-100", "3.1", "S/P1-C1.2"}},
    {{"100-619-300", "100-902-100", "3.2", "S/P1-C1.2"}},
    {{"100-619-300", "100-902-100", "3.3", "S/P1-C1.2"}},
    {{"100-619-300", "100-902-100", "4.0", "S/P1-C1.2"}},

    // 3.2 ("Goat"): S/P1-C1.3
    {{"100-619-300", "100-972-200", "3.2", "S/P1-C1.3"}},
    {{"100-619-300", "100-972-200", "3.3", "S/P1-C1.3"}},
    {{"100-619-300", "100-972-200", "4.0", "S/P1-C1.3"}},

    // 4.0 ("Seabiscuit"); S/P2-C2
    {{"100-862-200", "100-861-800", "4.0", "S/P2-C2"}}

};

} // namespace internal
} // namespace BAM
} // namespace PacBio

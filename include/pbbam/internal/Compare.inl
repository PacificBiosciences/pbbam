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
/// \file Compare.inl
/// \brief Inline implementations for the Compare class & inner classes.
//
// Author: Derek Barnett

#include "pbbam/Compare.h"

namespace PacBio {
namespace BAM {
namespace internal {

template <typename T, T> struct MemberFnProxy;

template<typename T, typename R, typename... Args, R (T::*fn)(Args...)const>
struct MemberFnProxy<R (T::*)(Args...)const, fn>
{
    static R call(const T& obj, Args&&... args)
    {
        return (obj.*fn)(std::forward<Args>(args)...);
    }
};

} // namespace internal

template<typename ValueType,
         typename Compare::MemberFunctionBaseHelper<ValueType>::MemberFnType fn,
         typename CompareType>
inline bool Compare::MemberFunctionBase<ValueType, fn, CompareType>::operator()(const BamRecord& lhs,
                                                                                const BamRecord& rhs) const
{
    using MemberFnType = typename Compare::MemberFunctionBaseHelper<ValueType>::MemberFnType;
    using Proxy = internal::MemberFnProxy<MemberFnType, fn>;

    CompareType cmp;
    return cmp(Proxy::call(lhs), Proxy::call(rhs));
}

inline bool Compare::None::operator()(const BamRecord&, const BamRecord&) const
{ return false; }

} // namespace BAM
} // namespace PacBio

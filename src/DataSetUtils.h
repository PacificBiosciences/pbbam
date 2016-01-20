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

#ifndef DATASETUTILS_H
#define DATASETUTILS_H

#include "pbbam/DataSetTypes.h"
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

namespace PacBio {
namespace BAM {
namespace internal {

static const std::string XML_VERSION = std::string { "3.0.1" };

template<typename T>
inline const T& NullObject(void)
{
    static const T empty;
    return empty;
}

template<>
inline const PacBio::BAM::DataSetMetadata& NullObject(void)
{
    static const PacBio::BAM::DataSetMetadata empty("", "");
    return empty;
}

inline
std::string GenerateUuid(void)
{
    static boost::uuids::random_generator gen;
    const boost::uuids::uuid uuid = gen();
    return boost::uuids::to_string(uuid);
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#ifndef FETCH_CHILD_CONST_REF
#define FETCH_CHILD_CONST_REF(Class, Type, Method) \
    \
    const PacBio::BAM::Type& Class::Method(void) const \
    { \
        try { \
            return Child<PacBio::BAM::Type>(#Type); \
        } catch (std::exception&) { \
            return internal::NullObject<PacBio::BAM::Type>(); \
        } \
    }
#endif

#ifndef FETCH_CHILD_REF
#define FETCH_CHILD_REF(Class, Type, Method) \
    \
    PacBio::BAM::Type& Class::Method(void) \
    { \
        if (!HasChild(#Type)) \
            AddChild(internal::NullObject<PacBio::BAM::Type>()); \
        return Child<PacBio::BAM::Type>(#Type); \
    }
#endif

#ifndef DEFINE_ACCESSORS
#define DEFINE_ACCESSORS(Class, Type, Method) \
    FETCH_CHILD_CONST_REF(Class, Type, Method) \
    FETCH_CHILD_REF(Class, Type, Method)
#endif

#endif // DATASETUTILS_H

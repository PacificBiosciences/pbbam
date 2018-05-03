// File Description
/// \file MD5.cpp
/// \brief Implements basic MD5 hash utilities
//
// Author: Brett Bowman

#include "PbbamInternalConfig.h"

#include "pbbam/MD5.h"

#include <stdexcept>

#include <htslib/hts.h>

namespace PacBio {
namespace BAM {

class Md5ContextHelper
{
public:
    Md5ContextHelper() : data_(hts_md5_init())
    {
        if (data_ == nullptr) throw std::runtime_error{"could not initialize MD5 context"};
    }

    ~Md5ContextHelper() { hts_md5_destroy(data_); }

public:
    std::string Encoded(const std::string& str)
    {
        hts_md5_update(data_, reinterpret_cast<void*>(const_cast<char*>(str.c_str())), str.size());

        unsigned char digest[16];
        hts_md5_final(digest, data_);

        char hexdigest[33];  // leave space for null-term
        hts_md5_hex(hexdigest, digest);

        return std::string{hexdigest, 32};
    }

private:
    hts_md5_context* data_;
};

/// \brief MD5 hash of a string as a 32-digit hexadecimal string
///
std::string MD5Hash(const std::string& str)
{
    Md5ContextHelper md5;
    return md5.Encoded(str);
}

}  // namespace BAM
}  // namespace PacBio

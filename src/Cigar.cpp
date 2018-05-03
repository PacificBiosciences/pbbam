// File Description
/// \file Cigar.cpp
/// \brief Implements the Cigar class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/Cigar.h"

#include <cstddef>
#include <cstdint>
#include <sstream>

namespace PacBio {
namespace BAM {

Cigar::Cigar(const std::string& cigarString) : std::vector<CigarOperation>{}
{
    size_t numberStart = 0;
    const size_t numChars = cigarString.size();
    for (size_t i = 0; i < numChars; ++i) {
        const char c = cigarString.at(i);
        if (!isdigit(c)) {
            const size_t distance = i - numberStart;
            const uint32_t length = stoul(cigarString.substr(numberStart, distance));
            push_back(CigarOperation(c, length));
            numberStart = i + 1;
        }
    }
}

std::string Cigar::ToStdString() const
{
    std::ostringstream s;
    const auto endIt = this->cend();
    for (auto iter = this->cbegin(); iter != endIt; ++iter) {
        const CigarOperation& cigar = (*iter);
        s << cigar.Length() << cigar.Char();
    }
    return s.str();
}

}  // namespace BAM
}  // namespace PacBio

// Author: Derek Barnett

#ifndef PBBAM_STRINGUTILITIES_H
#define PBBAM_STRINGUTILITIES_H

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

/// \brief Splits a string into tokens
///
/// \param[in] line     input string
/// \param[in] delim    character to split on
///
/// \returns vector of tokens
///
inline std::vector<std::string> Split(const std::string& line, const char delim = '\t')
{
    std::vector<std::string> tokens;
    std::istringstream lineStream(line);
    std::string token;
    while (std::getline(lineStream, token, delim))
        tokens.push_back(token);
    return tokens;
}

/// \brief Remove all whitespace from input string (start, end, & internal)
///
/// \param[in] input    original string
///
/// \returns new string with no whitespace
///
inline std::string RemoveAllWhitespace(std::string input)
{
    input.erase(
        std::remove_if(input.begin(), input.end(), [](const char c) { return std::isspace(c); }),
        input.end());
    return input;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_STRINGUTILITIES_H

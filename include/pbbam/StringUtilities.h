// Author: Derek Barnett

#ifndef PBBAM_STRINGUTILITIES_H
#define PBBAM_STRINGUTILITIES_H

#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

///
/// \brief Joins tokens into a single string
///
/// \param tokens   input strings
/// \param delim    delimiter character
///
/// \return joined string
///
std::string Join(const std::vector<std::string>& tokens, const char delim);

/// \brief Splits a string into tokens
///
/// \param[in] line     input string
/// \param[in] delim    character to split on
///
/// \returns vector of tokens
///
std::vector<std::string> Split(const std::string& line, const char delim = '\t');

/// \brief Remove all whitespace from input string (start, end, & internal)
///
/// \param[in] input    original string
///
/// \returns new string with no whitespace
///
std::string RemoveAllWhitespace(std::string input);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_STRINGUTILITIES_H

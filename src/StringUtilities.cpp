#include "PbbamInternalConfig.h"

#include <pbbam/StringUtilities.h>

#include <algorithm>
#include <sstream>

namespace PacBio {
namespace BAM {

std::string Join(const std::vector<std::string>& tokens, const char delim)
{
    std::string result;
    bool first = true;
    for (const auto& token : tokens) {
        if (!first) {
            result += delim;
        }
        result += token;
        first = false;
    }
    return result;
}

std::vector<std::string> Split(const std::string& line, const char delim)
{
    std::vector<std::string> tokens;
    std::istringstream lineStream(line);
    std::string token;
    while (std::getline(lineStream, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string RemoveAllWhitespace(std::string input)
{
    input.erase(
        std::remove_if(input.begin(), input.end(), [](const char c) { return std::isspace(c); }),
        input.end());
    return input;
}

}  // namespace BAM
}  // namespace PacBio

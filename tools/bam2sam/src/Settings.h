// Author: Derek Barnett

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <vector>

namespace bam2sam {

class Settings
{
public:
    Settings(void) : noHeader_(false), printHeaderOnly_(false) {}

public:
    std::string inputFilename_;
    bool noHeader_;
    bool printHeaderOnly_;
    std::vector<std::string> errors_;
};

}  // namespace bam2sam

#endif  // SETTINGS_H

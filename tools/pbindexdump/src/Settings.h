// Author: Derek Barnett

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <vector>

namespace pbindexdump {

class Settings
{
public:
    Settings(void) : format_("json"), jsonIndentLevel_(4), jsonRaw_(false) {}

public:
    std::string inputPbiFilename_;
    std::string format_;
    int jsonIndentLevel_;
    bool jsonRaw_;
    std::vector<std::string> errors_;
};

}  // namespace pbindexdump

#endif  // SETTINGS_H

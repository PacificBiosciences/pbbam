// Author: Derek Barnett

#ifndef PBINDEX_H
#define PBINDEX_H

#include <string>
#include <vector>

namespace pbindex {

class Settings
{
public:
    Settings();

public:
public:
    std::string inputBamFilename_;
    bool printPbiContents_;
    std::vector<std::string> errors_;
};

class PbIndex
{
public:
    static int Run(const Settings& settings);

private:
    static int Create(const Settings& settings);
    //    static int Print(const Settings& settings);
};

}  // namespace pbindex

#endif  // PBINDEX_H

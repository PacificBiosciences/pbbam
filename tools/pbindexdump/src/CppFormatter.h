// Author: Derek Barnett

#ifndef CPPFORMATTER_H
#define CPPFORMATTER_H

#include "IFormatter.h"

namespace pbindexdump {

class CppFormatter : public IFormatter
{
public:
    CppFormatter(const Settings& settings);
    void Run();
};

}  // namespace pbindexdump

#endif  // CPPFORMATTER_H

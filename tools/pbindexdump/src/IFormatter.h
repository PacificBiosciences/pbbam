// Author: Derek Barnett

#ifndef IFORMATTER_H
#define IFORMATTER_H

#include "Settings.h"

namespace pbindexdump {

class IFormatter
{
public:
    ~IFormatter(void) { }

public:
    virtual void Run(void) =0;

protected:
    const Settings& settings_;

protected:
    IFormatter(const Settings& settings)
        : settings_(settings)
    { }
};

} // namespace pbindexdump

#endif // IFORMATTER_H

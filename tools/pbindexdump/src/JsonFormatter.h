// Author: Derek Barnett

#ifndef JSONFORMATTER_H
#define JSONFORMATTER_H

#include "IFormatter.h"
#include "json.hpp"
#include <pbbam/PbiRawData.h>

namespace pbindexdump {

class JsonFormatter : public IFormatter
{
public:
    JsonFormatter(const Settings& settings);
    void Run();

private:
    void FormatMetadata();
    void FormatReferences();

    void FormatRaw();
    void FormatRecords();

    void Print();

private:
    PacBio::BAM::PbiRawData index_;
    nlohmann::json json_;
};

} // namespace pbindexdump

#endif // JSONFORMATTER_H

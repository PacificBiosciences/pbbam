// Author: Derek Barnett

#include "PbIndexDump.h"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

#include <pbbam/MakeUnique.h>

#include "CppFormatter.h"
#include "JsonFormatter.h"
using namespace pbindexdump;

void PbIndexDump::Run(const Settings& settings)
{
    std::unique_ptr<IFormatter> formatter(nullptr);
    if (settings.format_ == "json")
        formatter = std::make_unique<JsonFormatter>(settings);
    else if (settings.format_ == "cpp")
        formatter = std::make_unique<CppFormatter>(settings);
    else {
        std::string msg = {"unsupported output format requested: "};
        msg += settings.format_;
        throw std::runtime_error(msg);
    }
    assert(formatter);
    formatter->Run();
}

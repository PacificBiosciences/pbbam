// Author: Derek Barnett

#include "PbIndexDump.h"
#include <cassert>
#include "CppFormatter.h"
#include "JsonFormatter.h"
#include "pbbam/MakeUnique.h"
using namespace pbindexdump;
using namespace std;

void PbIndexDump::Run(const Settings& settings)
{
    std::unique_ptr<IFormatter> formatter(nullptr);
    if (settings.format_ == "json")
        formatter = std::make_unique<JsonFormatter>(settings);
    else if (settings.format_ == "cpp")
        formatter = std::make_unique<CppFormatter>(settings);
    else {
        string msg = {"unsupported output format requested: "};
        msg += settings.format_;
        throw runtime_error(msg);
    }
    assert(formatter);
    formatter->Run();
}

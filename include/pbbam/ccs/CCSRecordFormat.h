#ifndef PBBAM_CCS_CCSRECORDFORMAT_H
#define PBBAM_CCS_CCSRECORDFORMAT_H

#include "pbbam/Config.h"

#include <string>
#include <vector>

#include "pbbam/ccs/CCSHeader.h"
#include "pbbam/ccs/CCSRecord.h"

namespace PacBio {
namespace CCS {

struct CCSRecordFormat
{
    // header
    static CCSHeader DeserializeHeader(const std::vector<std::string>& lines);
    static std::vector<std::string> SerializeHeader(const CCSHeader& header);

    // record
    static CCSRecord DeserializeRecord(const std::string& line);
    static std::string SerializeRecord(const CCSRecord& record);
};

}  // namespace CCS
}  // namespace PacBio

#endif  //  PBBAM_CCS_CCSRECORDFORMAT_H

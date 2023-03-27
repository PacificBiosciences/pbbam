#ifndef PBBAM_CSV_CSVTYPES_H
#define PBBAM_CSV_CSVTYPES_H

#include <pbbam/Config.h>

#include <map>
#include <string>
#include <vector>

namespace PacBio {
namespace CSV {

using CsvHeader = std::vector<std::string>;

using CsvRecord = std::map<std::string, std::string>;

}  // namespace CSV
}  // namespace PacBio

#endif  // PBBAM_CSV_CSVTYPES_H

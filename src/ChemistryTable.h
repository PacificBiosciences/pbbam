#ifndef PBBAM_CHEMISTRYTABLE_H
#define PBBAM_CHEMISTRYTABLE_H

#include <pbbam/Config.h>

#include <array>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

using ChemistryTable = std::vector<std::array<std::string, 5>>;

const ChemistryTable& BuiltInChemistryTable();

const ChemistryTable& GetChemistryTableFromEnv();

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_CHEMISTRYTABLE_H

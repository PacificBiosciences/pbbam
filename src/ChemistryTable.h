// Author: Lance Hepler

#ifndef CHEMISTRYTABLE_H
#define CHEMISTRYTABLE_H

#include "pbbam/Config.h"

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

#endif  // CHEMISTRYTABLE_H

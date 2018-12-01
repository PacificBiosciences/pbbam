// Author: Lance Hepler

#ifndef CHEMISTRYTABLE_H
#define CHEMISTRYTABLE_H

#include <array>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

using ChemistryTable = std::vector<std::array<std::string, 4>>;

const ChemistryTable& BuiltInChemistryTable();

const ChemistryTable& GetChemistryTableFromEnv();

}  // namespace BAM
}  // namespace PacBio

#endif  // CHEMISTRYTABLE_H

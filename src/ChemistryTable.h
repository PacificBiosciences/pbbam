// Author: Lance Hepler

#ifndef CHEMISTRYTABLE_H
#define CHEMISTRYTABLE_H

#include <array>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {
namespace internal {

typedef std::vector<std::array<std::string, 4>> ChemistryTable;

extern const ChemistryTable BuiltInChemistryTable;

const ChemistryTable& GetChemistryTableFromEnv();

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // CHEMISTRYTABLE_H

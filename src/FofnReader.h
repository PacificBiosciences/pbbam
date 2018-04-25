// Author: Derek Barnett

#ifndef FOFNREADER_H
#define FOFNREADER_H

#include <iosfwd>
#include <string>
#include <vector>
#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {
namespace internal {

class FofnReader
{
public:
    static std::vector<std::string> Files(std::istream& in);
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // FOFNREADER_H

// Author: Derek Barnett

#ifndef FOFNREADER_H
#define FOFNREADER_H

#include <iosfwd>
#include <string>
#include <vector>
#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {

class FofnReader
{
public:
    static std::vector<std::string> Files(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FOFNREADER_H

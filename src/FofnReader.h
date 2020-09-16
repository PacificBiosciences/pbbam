#ifndef PBBAM_FOFNREADER_H
#define PBBAM_FOFNREADER_H

#include <pbbam/Config.h>

#include <iosfwd>
#include <string>
#include <vector>

#include <pbbam/DataSet.h>

namespace PacBio {
namespace BAM {

class FofnReader
{
public:
    static std::vector<std::string> Files(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FOFNREADER_H

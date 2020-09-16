#ifndef PBBAM_XMLREADER_H
#define PBBAM_XMLREADER_H

#include <pbbam/Config.h>

#include <iosfwd>
#include <memory>

#include <pbbam/DataSet.h>

namespace PacBio {
namespace BAM {

class XmlReader
{
public:
    static std::unique_ptr<DataSetBase> FromStream(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_XMLREADER_H

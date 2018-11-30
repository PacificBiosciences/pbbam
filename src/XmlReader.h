// Author: Derek Barnett

#ifndef XMLREADER_H
#define XMLREADER_H

#include <iosfwd>
#include <memory>

#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {

class XmlReader
{
public:
    static std::unique_ptr<DataSetBase> FromStream(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // XMLREADER_H

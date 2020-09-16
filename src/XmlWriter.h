#ifndef PBBAM_XMLWRITER_H
#define PBBAM_XMLWRITER_H

#include <pbbam/Config.h>

#include <iosfwd>
#include <memory>

#include <pbbam/DataSetTypes.h>

namespace PacBio {
namespace BAM {

class DataSetBase;

class XmlWriter
{
public:
    static void ToStream(const DataSetBase& dataset, std::ostream& out, DataSetPathMode pathMode);
    static void ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out,
                         DataSetPathMode pathMode);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_XMLWRITER_H

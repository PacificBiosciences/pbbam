#ifndef PBBAM_DATASETIO_H
#define PBBAM_DATASETIO_H

#include <pbbam/Config.h>

#include <pbbam/DataSet.h>

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class DataSetIO
{
public:
    // input
    static std::unique_ptr<DataSetBase> FromUri(const std::string& uri);
    static std::unique_ptr<DataSetBase> FromUris(const std::vector<std::string>& uris);
    static std::unique_ptr<DataSetBase> FromXmlString(const std::string& xml);

    // output
    static void ToFile(DataSetBase& dataset, const std::string& fn, DataSetPathMode pathMode);
    static void ToFile(const std::unique_ptr<DataSetBase>& dataset, const std::string& fn,
                       DataSetPathMode pathMode);
    static void ToStream(DataSetBase& dataset, std::ostream& out, DataSetPathMode pathMode);
    static void ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out,
                         DataSetPathMode pathMode);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_DATASETIO_H

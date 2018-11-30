// Author: Derek Barnett

#ifndef DATASETIO_H
#define DATASETIO_H

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
#include "pbbam/DataSet.h"

namespace PacBio {
namespace BAM {

class DataSetIO
{
public:
    // input
    static std::unique_ptr<DataSetBase> FromUri(const std::string& uri);
    static std::unique_ptr<DataSetBase> FromUris(const std::vector<std::string>& uris);

    static std::unique_ptr<DataSetBase> FromXmlString(const std::string& xml);

    //    static DataSetBase FromUri(const std::string& uri);
    //    static DataSetBase FromUris(const std::vector<std::string>& uris);

    //    // output
    static void ToFile(const std::unique_ptr<DataSetBase>& dataset, const std::string& fn);
    static void ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out);

    static void ToFile(const DataSetBase& dataset, const std::string& fn);
    static void ToStream(const DataSetBase& dataset, std::ostream& out);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // DATASETIO_H

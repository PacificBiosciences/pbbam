// Author: Derek Barnett

#ifndef XMLWRITER_H
#define XMLWRITER_H

#include <iosfwd>
#include <memory>

namespace PacBio {
namespace BAM {

class DataSetBase;

class XmlWriter
{
public:
    static void ToStream(const DataSetBase& dataset, std::ostream& out);
    static void ToStream(const std::unique_ptr<DataSetBase>& dataset, std::ostream& out);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // XMLWRITER_H

#ifndef PBBAM_CCS_CCSRECORDWRITER_H
#define PBBAM_CCS_CCSRECORDWRITER_H

#include "pbbam/Config.h"

#include <iosfwd>
#include <memory>
#include <string>

#include "pbbam/ccs/CCSHeader.h"
#include "pbbam/ccs/CCSRecord.h"

namespace PacBio {
namespace CCS {

///
/// Writes CCSRecords to stdout
///
class CCSRecordWriter
{
public:
    ///
    /// \brief Construct a new CCSRecordWriter object
    ///
    /// \param header
    ///
    CCSRecordWriter(const CCSHeader& header);

    CCSRecordWriter(const CCSHeader& header, std::ostream& out);

    ~CCSRecordWriter();

    ///
    /// \brief
    ///
    /// \param record
    /// \return true
    /// \return false
    ///
    void Write(const CCSRecord& record);

private:
    class CCSRecordWriterPrivate;
    std::unique_ptr<CCSRecordWriterPrivate> d_;
};

}  // namespace CCS
}  // namespace PacBio

#endif  //  PBBAM_CCS_CCSRECORDWRITER_H

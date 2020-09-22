#ifndef PBBAM_CCS_CCSRECORDREADER_H
#define PBBAM_CCS_CCSRECORDREADER_H

#include <pbbam/Config.h>

#include <iosfwd>
#include <memory>
#include <string>

#include <pbbam/ccs/CCSHeader.h>
#include <pbbam/ccs/CCSRecord.h>
#include <pbbam/internal/QueryBase.h>

namespace PacBio {
namespace CCS {

///
/// Reads CCSRecords from stdin
//
class CCSRecordReader : public BAM::internal::QueryBase<CCSRecord>
{
public:
    CCSRecordReader();
    CCSRecordReader(std::istream& in);
    ~CCSRecordReader();

public:
    const CCSHeader& Header() const;

    bool GetNext(CCSRecord& record);

private:
    class CCSRecordReaderPrivate;
    std::unique_ptr<CCSRecordReaderPrivate> d_;
};

}  // namespace CCS
}  // namespace PacBio

#endif  // PBBAM_CCS_CCSRECORD_H

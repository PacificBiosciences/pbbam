#include "PbbamInternalConfig.h"

#include <pbbam/ccs/CCSRecordWriter.h>

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

#include <pbbam/ccs/CCSRecordFormat.h>

namespace PacBio {
namespace CCS {

class CCSRecordWriter::CCSRecordWriterPrivate
{
public:
    CCSRecordWriterPrivate(const CCSHeader& header, std::ostream& out) : out_{out}
    {
        WriteHeader(header);
    }

    void WriteHeader(const CCSHeader& header)
    {
        const auto lines = CCSRecordFormat::SerializeHeader(header);
        for (const auto& line : lines)
            out_ << line << '\n';
        out_ << "#\n";
    }

    void Write(const CCSRecord& record)
    {
        out_ << CCSRecordFormat::SerializeRecord(record) << '\n';
    }

    std::ostream& out_;
};

CCSRecordWriter::CCSRecordWriter(const CCSHeader& header) : CCSRecordWriter{header, std::cout} {}

CCSRecordWriter::CCSRecordWriter(const CCSHeader& header, std::ostream& out)
    : d_{std::make_unique<CCSRecordWriterPrivate>(header, out)}
{
}

CCSRecordWriter::~CCSRecordWriter() = default;

void CCSRecordWriter::Write(const CCSRecord& record) { d_->Write(record); }

}  // namespace CCS
}  // namespace PacBio

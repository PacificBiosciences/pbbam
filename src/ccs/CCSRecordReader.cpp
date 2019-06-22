// File Description
/// \file CCSRecordReader.cpp
/// \brief Implements the CCSRecordReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/ccs/CCSRecordReader.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "pbbam/ccs/CCSRecordFormat.h"

namespace PacBio {
namespace CCS {

class CCSRecordReader::CCSRecordReaderPrivate
{
public:
    CCSRecordReaderPrivate(std::istream& in) : in_{in} { ReadHeader(); }

    void ReadHeader()
    {
        const std::string EndHeader{"#"};

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in_, line)) {
            if (line == EndHeader) break;
            lines.push_back(line);
        }
        header_ = CCSRecordFormat::DeserializeHeader(lines);
    }

    bool GetNext(CCSRecord& record)
    {
        if (!std::getline(in_, line_)) return false;  // indicates EOF
        record = CCSRecordFormat::DeserializeRecord(line_);
        return true;
    }

    std::istream& in_;
    std::string line_;
    CCSHeader header_;
};

CCSRecordReader::CCSRecordReader() : CCSRecordReader{std::cin} {}

CCSRecordReader::CCSRecordReader(std::istream& in)
    : d_{std::make_unique<CCSRecordReaderPrivate>(in)}
{
}

CCSRecordReader::~CCSRecordReader() = default;

bool CCSRecordReader::GetNext(CCSRecord& record) { return d_->GetNext(record); }

const CCSHeader& CCSRecordReader::Header() const { return d_->header_; }

}  // namespace CCS
}  // namespace PacBio
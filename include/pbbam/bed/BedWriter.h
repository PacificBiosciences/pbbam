// File Description
/// \file BedWriter.h
/// \brief Defines the BedWriter class.
//
// Author: Derek Barnett

#ifndef PBBAM_BEDWRITER_H
#define PBBAM_BEDWRITER_H

#include "pbbam/Config.h"

#include <memory>

namespace PacBio {
namespace BAM {

class GenomicInterval;

class BedWriter
{
public:
    explicit BedWriter(const std::string& fn);

    BedWriter(BedWriter&&) noexcept;
    BedWriter& operator=(BedWriter&&) noexcept;
    ~BedWriter();

public:
    void Write(const GenomicInterval& interval);

private:
    class BedWriterPrivate;
    std::unique_ptr<BedWriterPrivate> d_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // FASTQWRITER_H

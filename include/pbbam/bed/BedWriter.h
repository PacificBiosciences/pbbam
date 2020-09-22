#ifndef PBBAM_BED_BEDWRITER_H
#define PBBAM_BED_BEDWRITER_H

#include <pbbam/Config.h>

#include <memory>

#include <pbcopper/data/GenomicInterval.h>

namespace PacBio {
namespace BED {

class BedWriter
{
public:
    explicit BedWriter(const std::string& fn);

    BedWriter(BedWriter&&) noexcept;
    BedWriter& operator=(BedWriter&&) noexcept;
    ~BedWriter();

public:
    void Write(const Data::GenomicInterval& interval);

private:
    class BedWriterPrivate;
    std::unique_ptr<BedWriterPrivate> d_;
};

}  // namespace BED
}  // namespace PacBio

#endif  // PBBAM_BED_BEDWRITER_H

#include "PbbamInternalConfig.h"

#include <pbbam/bed/BedWriter.h>

#include <cassert>

#include <sstream>
#include <type_traits>

#include <pbbam/GenomicInterval.h>
#include <pbbam/TextFileWriter.h>

namespace PacBio {
namespace BED {

static_assert(!std::is_copy_constructible<BedWriter>::value,
              "BedWriter(const BedWriter&) is not = delete");
static_assert(!std::is_copy_assignable<BedWriter>::value,
              "BedWriter& operator=(const BedWriter&) is not = delete");

class BedWriter::BedWriterPrivate
{
public:
    explicit BedWriterPrivate(const std::string& filename) : writer_{filename} {}

    void Write(const Data::GenomicInterval& interval)
    {
        line_.str("");
        line_ << interval.Name() << '\t' << interval.Start() << '\t' << interval.Stop();
        writer_.Write(line_.str());
    }

private:
    std::ostringstream line_;
    BAM::TextFileWriter writer_;
};

BedWriter::BedWriter(const std::string& fn) : d_{std::make_unique<BedWriterPrivate>(fn)} {}

BedWriter::BedWriter(BedWriter&&) noexcept = default;

BedWriter& BedWriter::operator=(BedWriter&&) noexcept = default;

BedWriter::~BedWriter() = default;

void BedWriter::Write(const Data::GenomicInterval& interval) { d_->Write(interval); }

}  // namespace BED
}  // namespace PacBio

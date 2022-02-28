#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfWriter.h>

#include <pbbam/vcf/VcfFormat.h>
#include <pbbam/vcf/VcfHeader.h>
#include <pbbam/vcf/VcfVariant.h>
#include "../FileProducer.h"

#include <fstream>
#include <type_traits>

namespace PacBio {
namespace VCF {

struct VcfWriter::VcfWriterPrivate : public BAM::FileProducer
{
    VcfWriterPrivate(std::string fn, const VcfHeader& header)
        : BAM::FileProducer{std::move(fn)}, out_{TempFilename()}
    {
        out_ << VcfFormat::FormattedHeader(header) << '\n';
    }

    bool Write(const VcfVariant& var)
    {
        out_ << VcfFormat::FormattedVariant(var) << '\n';
        return true;  // TODO: handle errors
    }

    std::ofstream out_;
};

VcfWriter::VcfWriter(std::string fn, const VcfHeader& header)
    : d_{std::make_unique<VcfWriterPrivate>(std::move(fn), header)}
{}

VcfWriter::VcfWriter(VcfWriter&&) noexcept = default;

VcfWriter& VcfWriter::operator=(VcfWriter&&) noexcept = default;

VcfWriter::~VcfWriter() = default;

bool VcfWriter::Write(const VcfVariant& var) { return d_->Write(var); }

}  // namespace VCF
}  // namespace PacBio

// Author: Derek Barnett

#include <pbbam/vcf/VcfWriter.h>

#include <fstream>
#include <iostream>

#include <pbbam/MakeUnique.h>
#include <pbbam/vcf/VcfFormat.h>
#include <pbbam/vcf/VcfHeader.h>
#include <pbbam/vcf/VcfVariant.h>
#include "../FileProducer.h"

namespace PacBio {
namespace VCF {

struct VcfWriter::VcfWriterPrivate : public PacBio::BAM::internal::FileProducer
{
    VcfWriterPrivate(std::string fn, const VcfHeader& header)
        : PacBio::BAM::internal::FileProducer{std::move(fn)}, out_{TempFilename()}
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
{
}

bool VcfWriter::Write(const VcfVariant& var) { return d_->Write(var); }

VcfWriter::~VcfWriter() {}

}  // namespace VCF
}  // namespace PacBio

// Author: Derek Barnett

#include <pbbam/vcf/VcfReader.h>

namespace PacBio {
namespace VCF {

VcfReader::VcfReader(std::string fn) : VcfReader{VcfFile{std::move(fn)}} {}

VcfReader::VcfReader(const VcfFile& file) : in_{file.Filename()}, header_{file.Header()}
{
    // skip header lines
    const auto& header = file.Header();
    std::string line;
    for (size_t i = header.NumLines(); i > 0; --i)
        std::getline(in_, line);

    FetchNext();
}

void VcfReader::FetchNext()
{
    line_.clear();
    std::getline(in_, line_);
}

bool VcfReader::GetNext(VcfVariant& var)
{
    if (line_.empty()) return false;
    var = VcfVariant{line_};
    FetchNext();
    return true;
}

const VcfHeader& VcfReader::Header() const { return header_; }

}  // namespace VCF
}  // namespace PacBio

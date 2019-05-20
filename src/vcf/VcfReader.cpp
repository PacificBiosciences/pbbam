// Author: Derek Barnett

#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfReader.h>

namespace PacBio {
namespace VCF {

VcfReader::VcfReader(std::string fn) : VcfReader{VcfFile{std::move(fn)}} {}

VcfReader::VcfReader(const VcfFile& file)
    : in_{std::make_unique<std::ifstream>(file.Filename())}, header_{file.Header()}
{
    // skip header lines
    const auto& header = file.Header();
    std::string line;
    for (size_t i = header.NumLines(); i > 0; --i)
        std::getline(*in_, line);

    FetchNext();
}

VcfReader::VcfReader(VcfReader&&) noexcept = default;

VcfReader& VcfReader::operator=(VcfReader&&) PBBAM_NOEXCEPT_MOVE_ASSIGN = default;

VcfReader::~VcfReader() = default;

void VcfReader::FetchNext()
{
    line_.clear();
    std::getline(*in_, line_);
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

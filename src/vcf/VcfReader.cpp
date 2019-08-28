// Author: Derek Barnett

#include "../PbbamInternalConfig.h"

#include <pbbam/vcf/VcfReader.h>

#include <cassert>
#include <type_traits>

namespace PacBio {
namespace VCF {

static_assert(!std::is_copy_constructible<VcfReader>::value,
              "VcfReader(const VcfReader&) is not = delete");
static_assert(!std::is_copy_assignable<VcfReader>::value,
              "VcfReader& operator=(const VcfReader&) is not = delete");

static_assert(std::is_nothrow_move_constructible<VcfReader>::value ==
                  std::is_nothrow_move_constructible<std::ifstream>::value,
              "");
static_assert(std::is_nothrow_move_assignable<VcfReader>::value ==
                  std::is_nothrow_move_assignable<std::ifstream>::value,
              "");

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

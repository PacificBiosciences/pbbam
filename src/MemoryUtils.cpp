// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "MemoryUtils.h"

#include <cstdlib>
#include <cstring>
#include <string>

namespace PacBio {
namespace BAM {
namespace internal {

// -----------------
// BamHeaderMemory
// -----------------

BamHeader BamHeaderMemory::FromRawData(bam_hdr_t* hdr)
{
    // null input - error
    if (hdr == nullptr) throw std::runtime_error{"invalid BAM header"};

    // empty text input - ok
    if (hdr->text == nullptr || hdr->l_text == 0) return BamHeader();

    // parse normal SAM text input
    return BamHeader(std::string(hdr->text, hdr->l_text));
}

std::shared_ptr<bam_hdr_t> BamHeaderMemory::MakeRawHeader(const BamHeader& header)
{
    const std::string text = header.ToSam();
    std::shared_ptr<bam_hdr_t> rawData(sam_hdr_parse(text.size(), text.c_str()),
                                       internal::HtslibHeaderDeleter());
    rawData->ignore_sam_err = 0;
    rawData->cigar_tab = nullptr;
    rawData->l_text = text.size();
    rawData->text = static_cast<char*>(calloc(rawData->l_text + 1, 1));
    memcpy(rawData->text, text.c_str(), rawData->l_text);
    return rawData;
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

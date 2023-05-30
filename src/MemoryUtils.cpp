#include "PbbamInternalConfig.h"

#include "MemoryUtils.h"

#include <pbbam/Deleters.h>

#include <string>

#include <cstdlib>
#include <cstring>

namespace PacBio {
namespace BAM {

BamHeader BamHeaderMemory::FromRawData(bam_hdr_t* hdr)
{
    // null input - error
    if (hdr == nullptr) {
        throw std::runtime_error{"[pbbam] BAM header ERROR: null BAM header"};
    }

    // empty text input - ok
    if (hdr->text == nullptr || hdr->l_text == 0) {
        return BamHeader();
    }

    // parse normal SAM text input
    return BamHeader(std::string(hdr->text, hdr->l_text));
}

std::shared_ptr<bam_hdr_t> BamHeaderMemory::MakeRawHeader(const BamHeader& header)
{
    const std::string text = header.ToSam();
    std::shared_ptr<bam_hdr_t> rawData(sam_hdr_parse(text.size(), text.c_str()),
                                       HtslibHeaderDeleter());
    rawData->ignore_sam_err = 0;
    rawData->l_text = text.size();
    rawData->text = static_cast<char*>(std::calloc(rawData->l_text + 1, 1));
    std::memcpy(rawData->text, text.c_str(), rawData->l_text);

// HTS_VERSION only added >= v1.10, and this step is only necessary before then
#ifndef HTS_VERSION
    rawData->cigar_tab = nullptr;
#endif

    return rawData;
}

}  // namespace BAM
}  // namespace PacBio

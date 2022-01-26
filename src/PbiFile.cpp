#include "PbbamInternalConfig.h"

#include <pbbam/PbiFile.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/PbiBuilder.h>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

void PbiFile::CreateFrom(const BamFile& bamFile,
                         const PbiBuilder::CompressionLevel compressionLevel,
                         const size_t numThreads)
{
    PbiBuilder builder{bamFile.PacBioIndexFilename(), bamFile.Header().Sequences().size(),
                       compressionLevel, numThreads};
    BamReader reader{bamFile};
    BamRecord b;
    int64_t offset = reader.VirtualTell();
    while (reader.GetNext(b)) {
        builder.AddRecord(b, offset);
        offset = reader.VirtualTell();
    }
}

}  // namespace BAM
}  // namespace PacBio

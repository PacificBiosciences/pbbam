#ifndef PBBAM_PBIFILE_H
#define PBBAM_PBIFILE_H

#include <pbbam/Config.h>

#include <pbbam/PbiBuilder.h>

#include <string>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

class BamFile;

struct PbiFile
{

    /// \brief This enum describes the PBI file sections
    ///
    enum Section
    {
        BASIC = 0x0000,      ///< BasicData     (required)
        MAPPED = 0x0001,     ///< MappedData    (always optional)
        REFERENCE = 0x0002,  ///< ReferenceData (always optional)
        BARCODE = 0x0004,    ///< BarcodeData   (always optional)

        ALL = BASIC | MAPPED | REFERENCE | BARCODE  ///< Synonym for 'all sections'
    };

    /// \brief Helper typedef for storing multiple Section flags.
    ///
    using Sections = uint16_t;

    /// \brief This enum describes the PBI file version.
    enum VersionEnum
    {
        Version_3_0_0 = 0x030000,  ///< v3.0.0
        Version_3_0_1 = 0x030001,  ///< v3.0.1
        Version_3_0_2 = 0x030002,  ///< v3.0.2
        Version_4_0_0 = 0x040000,  ///< v4.0.0

        CurrentVersion = Version_4_0_0  ///< Synonym for the current PBI version.
    };

    ///
    /// \brief The BasicField enum
    ///
    enum class BasicField
    {
        RG_ID,
        Q_START,
        Q_END,
        ZMW,
        READ_QUALITY,
        CONTEXT_FLAG,
        VIRTUAL_OFFSET
    };

    ///
    /// \brief The MappedField enum
    ///
    enum class MappedField
    {
        T_ID,
        T_START,
        T_END,
        A_START,
        A_END,
        N_M,
        N_MM,
        N_INS,
        N_DEL,
        MAP_QUALITY,
        STRAND,
        N_INS_OPS,
        N_DEL_OPS,
    };

    ///
    /// \brief The BarcodeField enum
    ///
    enum class BarcodeField
    {
        BC_FORWARD,
        BC_REVERSE,
        BC_QUALITY
    };

    /// \brief Builds PBI index data from the supplied %BAM file and writes a
    ///        ".pbi" file.
    ///
    /// \param[in] bamFile source %BAM file
    ///
    /// \throws std::runtime_error if index file could not be created
    ///
    static void CreateFrom(
        const BamFile& bamFile,
        PbiBuilder::CompressionLevel compressionLevel = PbiBuilder::DefaultCompression,
        std::size_t numThreads = 4);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBIFILE_H

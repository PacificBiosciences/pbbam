#ifndef PBBAM_RUNMETADATAPARSER_H
#define PBBAM_RUNMETADATAPARSER_H

#include <pbbam/Config.h>

#include <pbbam/RunMetadata.h>
#include "pugixml/pugixml.hpp"

#include <iosfwd>
#include <map>
#include <string>

namespace PacBio {
namespace BAM {

// clang-format off
namespace Element {

constexpr char AUTOMATION[]              = "Automation";
constexpr char AUTOMATION_PARAMETER[]    = "AutomationParameter";
constexpr char AUTOMATION_PARAMETERS[]   = "AutomationParameters";
constexpr char BINDING_KIT[]             = "BindingKit";
constexpr char CELL_NFC_INDEX[]          = "CellNFCIndex";
constexpr char COLLECTION_METADATA[]     = "CollectionMetadata";
constexpr char COLLECTION_NUMBER[]       = "CollectionNumber";
constexpr char COLLECTIONS[]             = "Collections";
constexpr char CONTROL_KIT[]             = "ControlKit";
constexpr char CUSTOM_SEQUENCE[]         = "CustomSequence";
constexpr char DATASET_METADATA[]        = "DataSetMetadata";
constexpr char EXPERIMENT_CONTAINER[]    = "ExperimentContainer";
constexpr char EXPOSURE[]                = "Exposure";
constexpr char EXTEND_FIRST[]            = "ExtendFirst";
constexpr char EXTENSION_TIME[]          = "ExtensionTime";
constexpr char EXTRA_IM_WASHES[]         = "ExtraIMWashes";
constexpr char HAS_N2_SWITCH[]           = "HasN2Switch";
constexpr char HQRF_METHOD[]             = "HQRFMethod";
constexpr char IMMOBILIZATION_TIME[]     = "ImmobilizationTime";
constexpr char INSERT_SIZE[]             = "InsertSize";
constexpr char LEFT_ADAPTER[]            = "LeftAdapter";
constexpr char LEFT_ADAPTOR_SEQUENCE[]   = "LeftAdaptorSequence";
constexpr char LEFT_PRIMER_SEQUENCE[]    = "LeftPrimerSequence";
constexpr char MOVIE_LENGTH[]            = "MovieLength";
constexpr char OUTPUTS[]                 = "Outputs";
constexpr char PACBIO_DATA_MODEL[]       = "PacBioDataModel";
constexpr char PART_NUMBER[]             = "PartNumber";
constexpr char PCD_IN_PLATE[]            = "PCDinPlate";
constexpr char PRE_EXTENSION_WORKFLOW[]  = "PreExtensionWorkflow";
constexpr char RIGHT_ADAPTER[]           = "RightAdapter";
constexpr char RIGHT_ADAPTOR_SEQUENCE[]  = "RightAdaptorSequence";
constexpr char RIGHT_PRIMER_SEQUENCE[]   = "RightPrimerSequence";
constexpr char RUN[]                     = "Run";
constexpr char RUNS[]                    = "Runs";
constexpr char SEQUENCE[]                = "Sequence";
constexpr char SEQUENCING_KIT_PLATE[]    = "SequencingKitPlate";
constexpr char SNR_CUT[]                 = "SNRCut";
constexpr char SUBREADSET[]              = "SubreadSet";
constexpr char SUBREADSETS[]             = "SubreadSets";
constexpr char TEMPLATE_PREP_KIT[]       = "TemplatePrepKit";
constexpr char TIP_SEARCH_MAX_DURATION[] = "TipSearchMaxDuration";
constexpr char USE_STAGE_HOT_START[]     = "UseStageHotStart";

}  // namespace Element
// clang-format on

struct RunMetadataParser
{
    static CollectionMetadata LoadCollection(const std::string& metadataXmlFn);
    static CollectionMetadata LoadCollection(std::istream& in);

    static std::map<std::string, CollectionMetadata> LoadCollections(
        const std::string& runMetadataXmlFn);
    static std::map<std::string, CollectionMetadata> LoadCollections(std::istream& in);
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_RUNMETADATAPARSER_H

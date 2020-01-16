// Author: Derek Barnett

#ifndef RUNMETADATAPARSER_H
#define RUNMETADATAPARSER_H

#include "pbbam/Config.h"

#include <iosfwd>
#include <map>
#include <string>

#include "pbbam/RunMetadata.h"

#include "pugixml/pugixml.hpp"

namespace PacBio {
namespace BAM {

// clang-format off
namespace Element {

constexpr const char Automation[]           = "Automation";
constexpr const char AutomationParameter[]  = "AutomationParameter";
constexpr const char AutomationParameters[] = "AutomationParameters";
constexpr const char BindingKit[]           = "BindingKit";
constexpr const char CellNFCIndex[]         = "CellNFCIndex";
constexpr const char CollectionMetadata[]   = "CollectionMetadata";
constexpr const char CollectionNumber[]     = "CollectionNumber";
constexpr const char Collections[]          = "Collections";
constexpr const char ControlKit[]           = "ControlKit";
constexpr const char CustomSequence[]       = "CustomSequence";
constexpr const char DataSetMetadata[]      = "DataSetMetadata";
constexpr const char ExperimentContainer[]  = "ExperimentContainer";
constexpr const char Exposure[]             = "Exposure";
constexpr const char ExtendFirst[]          = "ExtendFirst";
constexpr const char ExtensionTime[]        = "ExtensionTime";
constexpr const char ExtraIMWashes[]        = "ExtraIMWashes";
constexpr const char HasN2Switch[]          = "HasN2Switch";
constexpr const char HQRFMethod[]           = "HQRFMethod";
constexpr const char ImmobilizationTime[]   = "ImmobilizationTime";
constexpr const char InsertSize[]           = "InsertSize";
constexpr const char LeftAdapter[]          = "LeftAdapter";
constexpr const char LeftAdaptorSequence[]  = "LeftAdaptorSequence";
constexpr const char LeftPrimerSequence[]   = "LeftPrimerSequence";
constexpr const char MovieLength[]          = "MovieLength";
constexpr const char Outputs[]              = "Outputs";
constexpr const char PacBioDataModel[]      = "PacBioDataModel";
constexpr const char PartNumber[]           = "PartNumber";
constexpr const char PCDinPlate[]           = "PCDinPlate";
constexpr const char PreExtensionWorkflow[] = "PreExtensionWorkflow";
constexpr const char RightAdapter[]         = "RightAdapter";
constexpr const char RightAdaptorSequence[] = "RightAdaptorSequence";
constexpr const char RightPrimerSequence[]  = "RightPrimerSequence";
constexpr const char Run[]                  = "Run";
constexpr const char Runs[]                 = "Runs";
constexpr const char Sequence[]             = "Sequence";
constexpr const char SequencingKitPlate[]   = "SequencingKitPlate";
constexpr const char SNRCut[]               = "SNRCut";
constexpr const char SubreadSet[]           = "SubreadSet";
constexpr const char SubreadSets[]          = "SubreadSets";
constexpr const char TemplatePrepKit[]      = "TemplatePrepKit";
constexpr const char TipSearchMaxDuration[] = "TipSearchMaxDuration";
constexpr const char UseStageHotStart[]     = "UseStageHotStart";

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

#endif  // RUNMETADATAPARSER_H
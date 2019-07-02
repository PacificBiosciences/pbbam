// Author: Lance Hepler

#include "PbbamInternalConfig.h"

#include "ChemistryTable.h"

#include <cstdlib>
#include <fstream>
#include <map>

#include "FileUtils.h"
#include "pbbam/exception/BundleChemistryMappingException.h"
#include "pugixml/pugixml.hpp"

namespace PacBio {
namespace BAM {
namespace {

ChemistryTable ChemistryTableFromXml(const std::string& mappingXml)
{
    if (!FileUtils::Exists(mappingXml))
        throw BundleChemistryMappingException{
            mappingXml, "SMRT_CHEMISTRY_BUNDLE_DIR defined but file not found"};

    std::ifstream in(mappingXml);
    pugi::xml_document doc;
    const pugi::xml_parse_result loadResult = doc.load(in);
    if (loadResult.status != pugi::status_ok)
        throw BundleChemistryMappingException{
            mappingXml, "unparseable XML, error code:" + std::to_string(loadResult.status)};

    // parse top-level attributes
    pugi::xml_node rootNode = doc.document_element();
    if (rootNode == pugi::xml_node())
        throw BundleChemistryMappingException{mappingXml, "could not fetch XML root node"};

    if (std::string(rootNode.name()) != "MappingTable")
        throw BundleChemistryMappingException{mappingXml, "MappingTable not found"};

    ChemistryTable table;
    try {
        for (const auto& childNode : rootNode) {
            const std::string childName = childNode.name();
            if (childName != "Mapping") continue;
            table.emplace_back(
                std::array<std::string, 4>{{childNode.child("BindingKit").child_value(),
                                            childNode.child("SequencingKit").child_value(),
                                            childNode.child("SoftwareVersion").child_value(),
                                            childNode.child("SequencingChemistry").child_value()}});
        }
    } catch (std::exception& e) {
        const std::string msg = std::string{"Mapping entries unparseable - "} + e.what();
        throw BundleChemistryMappingException{mappingXml, msg};
    }
    return table;
}

}  // namespace

const ChemistryTable& BuiltInChemistryTable()
{
    // clang-format off
    static ChemistryTable builtin{

        // BindingKit, SequencingKit, BasecallerVersion, Chemistry

        // 3.0 ("Dromedary"): S/P1-C1/beta
        {{"100-619-300", "100-620-000", "3.0", "S/P1-C1/beta"}},
        {{"100-619-300", "100-620-000", "3.1", "S/P1-C1/beta"}},

        // 3.1 ("Echidna"): S/P1-C1.1
        {{"100-619-300", "100-867-300", "3.1", "S/P1-C1.1"}},
        {{"100-619-300", "100-867-300", "3.2", "S/P1-C1.1"}},
        {{"100-619-300", "100-867-300", "3.3", "S/P1-C1.1"}},

        // 3.1.1 ("Flea"): S/P1-C1.2
        {{"100-619-300", "100-902-100", "3.1", "S/P1-C1.2"}},
        {{"100-619-300", "100-902-100", "3.2", "S/P1-C1.2"}},
        {{"100-619-300", "100-902-100", "3.3", "S/P1-C1.2"}},
        {{"100-619-300", "100-902-100", "4.0", "S/P1-C1.2"}},
        {{"100-619-300", "100-902-100", "4.1", "S/P1-C1.2"}},

        // 3.2 ("Goat"): S/P1-C1.3
        {{"100-619-300", "100-972-200", "3.2", "S/P1-C1.3"}},
        {{"100-619-300", "100-972-200", "3.3", "S/P1-C1.3"}},
        {{"100-619-300", "100-972-200", "4.0", "S/P1-C1.3"}},
        {{"100-619-300", "100-972-200", "4.1", "S/P1-C1.3"}},

        // 4.0 ("Seabiscuit"); S/P2-C2
        {{"100-862-200", "100-861-800", "4.0", "S/P2-C2"}},
        {{"100-862-200", "100-861-800", "4.1", "S/P2-C2"}},
        {{"100-862-200", "101-093-700", "4.1", "S/P2-C2"}},

        // 5.0 ("Iguana"); S/P2-C2
        {{"100-862-200", "100-861-800", "5.0", "S/P2-C2/5.0"}},
        {{"100-862-200", "101-093-700", "5.0", "S/P2-C2/5.0"}},

        // 5.0.1 ChemRel ("Sequel® Sequencing Plate Silwet"); S/P2-C2
        {{"100-862-200", "101-309-500", "5.0", "S/P2-C2/5.0"}},
        // 5.0.1 ChemRel ("Sequel® Sequencing Plate Silwet (4 rxn)"); S/P2-C2
        {{"100-862-200", "101-309-400", "5.0", "S/P2-C2/5.0"}},

        // --- SG1/16509P/PA5.0 ---
        // 2.1 binding kit/5.1PA support with ..
        // 5.0 ("Iguana"); S/P2-C2
        {{"101-365-900", "100-861-800", "5.0", "S/P2-C2/5.0"}},
        {{"101-365-900", "101-093-700", "5.0", "S/P2-C2/5.0"}},

        // 5.0.1 ChemRel; Sequel® Binding Kit 2.1; S/P2-C2
        {{"101-365-900", "101-309-500", "5.0", "S/P2-C2/5.0"}}, // Sequel® Sequencing Plate 2.1 Silwet (8 rxn)
        {{"101-365-900", "101-309-400", "5.0", "S/P2-C2/5.0"}}, // Sequel® Sequencing Plate 2.1 Silwet (4 rxn)

        // 5.0.1 ChemRel; Sequel® Binding Kit 3.0; S/P3-C3
        {{"101-500-400", "101-427-500", "5.0", "S/P3-C3/5.0"}}, // Sequel® Sequencing Plate 3.0 (8 rxn)
        {{"101-500-400", "101-427-800", "5.0", "S/P3-C3/5.0"}}, // Sequel® Sequencing Plate 3.0 (4 rxn)

        // 5.0.1 ChemRel; Sequel® Dev Binding Kit; S/P2-C2
        {{"101-490-800", "101-490-900", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel II® Sequencing Plate (4 rxn)
        {{"101-490-800", "101-491-000", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel II® Sequencing Plate (8 rxn)

        // 5.0.1 ChemRel; Sequel® Sequencing Plate 3.1 for Dynamic Loading placeholder (4 rxn)
        {{"101-500-400", "101-646-800", "5.0", "S/P3-C3/5.0"}}, // Sequel® Sequencing Plate 3.1 for Dynamic Loading placeholder

        // 5.0.1 ChemRel; Sequel® Dev Sequencing Plate Dyn Loading (4 rxn)
        {{"101-490-800", "101-644-500", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel® Dev Sequencing Plate Dyn Loading

        // 5.0.1 ChemRel; Sequel® Sequencing Plate Dyn Loading (4 rxn)
        {{"101-490-800", "101-717-100", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel® Dev Sequencing Plate Dyn Loading

        // 5.0.1 ChemRel; Sequel® Dev Sequencing Plate Dyn Loading (4 rxn)
        {{"101-717-300", "101-644-500", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel® Dev Sequencing Plate Dyn Loading
        // 5.0.1 ChemRel; Sequel® Sequencing Plate Dyn Loading (4 rxn)
        {{"101-717-300", "101-717-100", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel® Dev Sequencing Plate Dyn Loading

        // 5.0.1 ChemRel; Sequel® Dev Sequencing Plate Dyn Loading (4 rxn)
        {{"101-717-400", "101-644-500", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel® Dev Sequencing Plate Dyn Loading
        // 5.0.1 ChemRel; Sequel® Sequencing Plate Dyn Loading (4 rxn)
        {{"101-717-400", "101-717-100", "5.0", "S/P3-C1/5.0-8M"}}, // Sequel® Dev Sequencing Plate Dyn Loading

        // Sequel® II Binding Kit 2.0EA; Sequel® II Sequencing Plate 2.0EA (4 Rxn)
        {{"101-789-500", "101-789-300", "5.0", "S/P4-C2/5.0-8M"}}
    };
    // clang-format on

    return builtin;
}

const ChemistryTable& GetChemistryTableFromEnv()
{
    static const ChemistryTable empty{};
    static std::map<std::string, ChemistryTable> tableCache;

    std::string chemPath;
    const char* pth = getenv("SMRT_CHEMISTRY_BUNDLE_DIR");
    if (pth != nullptr && pth[0] != '\0')
        chemPath = pth;
    else
        return empty;

    auto it = tableCache.find(chemPath);
    if (it != tableCache.end()) return it->second;

    auto tbl = ChemistryTableFromXml(chemPath + "/chemistry.xml");
    it = tableCache.emplace(std::move(chemPath), std::move(tbl)).first;
    return it->second;
}

}  // namespace BAM
}  // namespace PacBio

// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Lance Hepler

#include "pbbam/exception/BundleChemistryMappingException.h"
#include "ChemistryTable.h"
#include "FileUtils.h"
#include "pugixml/pugixml.hpp"
#include <fstream>
#include <map>
#include <cstdlib>

namespace PacBio {
namespace BAM {
namespace internal {

extern const ChemistryTable BuiltInChemistryTable = {

    // BindingKit, SequencingKit, BasecallerVersion, Chemistry

    // RS
    {{"100356300",   "100356200",   "2.1", "P6-C4"}},
    {{"100356300",   "100356200",   "2.3", "P6-C4"}},
    {{"100356300",   "100612400",   "2.1", "P6-C4"}},
    {{"100356300",   "100612400",   "2.3", "P6-C4"}},
    {{"100372700",   "100356200",   "2.1", "P6-C4"}},
    {{"100372700",   "100356200",   "2.3", "P6-C4"}},
    {{"100372700",   "100612400",   "2.1", "P6-C4"}},
    {{"100372700",   "100612400",   "2.3", "P6-C4"}},

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
    {{"100-862-200", "101-093-700", "5.0", "S/P2-C2/5.0"}}

};

ChemistryTable ChemistryTableFromXml(const std::string& mappingXml)
{
    if (!FileUtils::Exists(mappingXml))
        throw BundleChemistryMappingException(
                mappingXml, "SMRT_CHEMISTRY_BUNDLE_DIR defined but file not found");

    std::ifstream in(mappingXml);
    pugi::xml_document doc;
    const pugi::xml_parse_result& loadResult = doc.load(in);
    if (loadResult.status != pugi::status_ok)
        throw BundleChemistryMappingException(
                mappingXml, std::string("unparseable XML, error code:") + std::to_string(loadResult.status));

    // parse top-level attributes
    pugi::xml_node rootNode = doc.document_element();
    if (rootNode == pugi::xml_node())
        throw BundleChemistryMappingException(mappingXml, "could not fetch XML root node");

    if (std::string(rootNode.name()) != "MappingTable")
        throw BundleChemistryMappingException(mappingXml, "MappingTable not found");

    ChemistryTable table;
    try {
        for (const auto& childNode : rootNode) {
            const std::string childName = childNode.name();
            if (childName != "Mapping") continue;
            table.emplace_back(std::array<std::string, 4>{{
                childNode.child("BindingKit").child_value(),
                childNode.child("SequencingKit").child_value(),
                childNode.child("SoftwareVersion").child_value(),
                childNode.child("SequencingChemistry").child_value()}});
        }
    } catch(std::exception& e) {
        const std::string msg = std::string{"Mapping entries unparseable - "} + e.what();
        throw BundleChemistryMappingException(mappingXml, msg);
    }
    return table;
}

const ChemistryTable& GetChemistryTableFromEnv()
{
    static const ChemistryTable empty{};
    static std::map<std::string, ChemistryTable> tableCache;

    std::string chemPath;
    const char* pth = getenv("SMRT_CHEMISTRY_BUNDLE_DIR");
    if (pth != nullptr && pth[0] != '\0')
        chemPath = pth;
    else return empty;

    auto it = tableCache.find(chemPath);
    if (it != tableCache.end()) return it->second;

    auto tbl = ChemistryTableFromXml(chemPath + "/chemistry.xml");
    it = tableCache.emplace(std::move(chemPath), std::move(tbl)).first;
    return it->second;
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

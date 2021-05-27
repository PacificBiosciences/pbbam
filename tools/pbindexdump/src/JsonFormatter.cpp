#include "JsonFormatter.h"

#include <cmath>
#include <cstdint>

#include <iostream>
#include <stdexcept>
#include <string>

#include <pbbam/PbiFile.h>
#include <pbbam/PbiRawData.h>

#include <pbcopper/json/JSON.h>

namespace PacBio {
namespace PbIndexDump {
namespace {

void FormatMetadata(const BAM::PbiRawData& index, JSON::Json& result)
{
    std::string version;
    switch (index.Version()) {
        case BAM::PbiFile::Version_3_0_0:
            version = "3.0.0";
            break;
        case BAM::PbiFile::Version_3_0_1:
            version = "3.0.1";
            break;
        case BAM::PbiFile::Version_3_0_2:
            version = "3.0.2";
            break;
        case BAM::PbiFile::Version_4_0_0:
            version = "4.0.0";
            break;
        default:
            throw std::runtime_error{"unsupported PBI version encountered"};
    }

    JSON::Json fileSections;
    fileSections.push_back("BasicData");
    if (index.HasBarcodeData()) {
        fileSections.push_back("BarcodeData");
    }
    if (index.HasMappedData()) {
        fileSections.push_back("MappedData");
    }
    if (index.HasReferenceData()) {
        fileSections.push_back("ReferenceData");
    }

    result["version"] = version;
    result["fileSections"] = fileSections;
    result["numReads"] = index.NumReads();
}

void FormatRaw(const BAM::PbiRawData& index, JSON::Json& result)
{
    const BAM::PbiRawBasicData& basicData = index.BasicData();
    result["basicData"]["rgId"] = basicData.rgId_;
    result["basicData"]["qStart"] = basicData.qStart_;
    result["basicData"]["qEnd"] = basicData.qEnd_;
    result["basicData"]["holeNumber"] = basicData.holeNumber_;
    result["basicData"]["readQual"] = basicData.readQual_;
    result["basicData"]["ctxtFlag"] = basicData.ctxtFlag_;
    result["basicData"]["fileOffset"] = basicData.fileOffset_;

    if (index.HasBarcodeData()) {
        const BAM::PbiRawBarcodeData& barcodeData = index.BarcodeData();
        result["barcodeData"]["bcForward"] = barcodeData.bcForward_;
        result["barcodeData"]["bcReverse"] = barcodeData.bcReverse_;
        result["barcodeData"]["bcQuality"] = barcodeData.bcQual_;
    }

    if (index.HasMappedData()) {
        const BAM::PbiRawMappedData& mappedData = index.MappedData();

        // casts to force -1 if unmapped
        result["mappedData"]["tId"] = mappedData.tId_;
        result["mappedData"]["tStart"] = mappedData.tStart_;
        result["mappedData"]["tEnd"] = mappedData.tEnd_;

        result["mappedData"]["aStart"] = mappedData.aStart_;
        result["mappedData"]["aEnd"] = mappedData.aEnd_;
        result["mappedData"]["revStrand"] = mappedData.revStrand_;
        result["mappedData"]["nM"] = mappedData.nM_;
        result["mappedData"]["nMM"] = mappedData.nMM_;
        result["mappedData"]["mapQV"] = mappedData.mapQV_;

        if (mappedData.hasIndelOps_) {
            result["mappedData"]["nInsOps"] = mappedData.nInsOps_;
            result["mappedData"]["nDelOps"] = mappedData.nDelOps_;
        }
    }
}

void FormatRecords(const BAM::PbiRawData& index, JSON::Json& result)
{
    JSON::Json reads;
    const uint32_t numReads = index.NumReads();
    const bool hasBarcodeData = index.HasBarcodeData();
    const bool hasMappedData = index.HasMappedData();

    for (uint32_t i = 0; i < numReads; ++i) {

        JSON::Json read;

        // common data
        const BAM::PbiRawBasicData& basicData = index.BasicData();
        read["rgId"] = basicData.rgId_[i];
        read["qStart"] = basicData.qStart_[i];
        read["qEnd"] = basicData.qEnd_[i];
        read["holeNumber"] = basicData.holeNumber_[i];
        read["readQuality"] = basicData.readQual_[i];
        read["contextFlag"] = basicData.ctxtFlag_[i];
        read["fileOffset"] = basicData.fileOffset_[i];

        // barcode data, if present
        if (hasBarcodeData) {
            const BAM::PbiRawBarcodeData& barcodeData = index.BarcodeData();
            read["bcForward"] = barcodeData.bcForward_[i];
            read["bcReverse"] = barcodeData.bcReverse_[i];
            read["bcQuality"] = barcodeData.bcQual_[i];
        }

        // mapping data, if present
        if (hasMappedData) {
            const BAM::PbiRawMappedData& mappedData = index.MappedData();

            // casts to force -1 if unmapped
            read["tId"] = static_cast<int32_t>(mappedData.tId_[i]);
            read["tStart"] = static_cast<int32_t>(mappedData.tStart_[i]);
            read["tEnd"] = static_cast<int32_t>(mappedData.tEnd_[i]);

            read["aStart"] = mappedData.aStart_[i];
            read["aEnd"] = mappedData.aEnd_[i];
            read["nM"] = mappedData.nM_[i];
            read["nMM"] = mappedData.nMM_[i];
            read["mapQuality"] = mappedData.mapQV_[i];
            read["reverseStrand"] = mappedData.revStrand_[i];

            if (mappedData.hasIndelOps_) {
                read["nInsOps"] = mappedData.nInsOps_[i];
                read["nDelOps"] = mappedData.nDelOps_[i];
            }
        }

        reads.push_back(std::move(read));
    }
    result["reads"] = reads;
}

void FormatReferences(const BAM::PbiRawData& index, JSON::Json& result)
{
    if (index.HasReferenceData()) {
        JSON::Json references;
        const auto& referenceData = index.ReferenceData();
        for (const auto& entry : referenceData.entries_) {
            JSON::Json element;
            element["tId"] = static_cast<int32_t>(entry.tId_);
            element["beginRow"] = static_cast<int32_t>(entry.beginRow_);
            element["endRow"] = static_cast<int32_t>(entry.endRow_);
            references.push_back(std::move(element));
        }
        result["references"] = references;
    }
}
}  // namespace

void JsonFormatter::Run(const Settings& settings)
{
    const BAM::PbiRawData index{settings.InputFile};
    JSON::Json result;

    FormatMetadata(index, result);
    FormatReferences(index, result);

    if (settings.JsonRaw) {
        FormatRaw(index, result);
    } else {
        FormatRecords(index, result);
    }

    // print
    std::cout << result.dump(settings.JsonIndentLevel) << '\n';
}

}  // namespace PbIndexDump
}  // namespace PacBio

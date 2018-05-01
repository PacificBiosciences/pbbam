// Author: Derek Barnett

#include "CppFormatter.h"
#include <pbbam/PbiRawData.h>

#include <cstdint>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace pbindexdump;

namespace pbindexdump {

static std::string printCppReferenceData(const PacBio::BAM::PbiRawReferenceData& referenceData)
{
    auto result = std::string{""};
    for (const PacBio::BAM::PbiReferenceEntry& entry : referenceData.entries_) {
        if (!result.empty()) result.append(",\n");
        result.append(std::string{"    PbiReferenceEntry{"} + std::to_string(entry.tId_) + "," +
                      std::to_string(entry.beginRow_) + "," + std::to_string(entry.endRow_) +
                      std::string{"}"});
    }
    if (!result.empty()) result.append("\n");
    return result;
}

template <typename T>
std::string printVectorElements(const std::vector<T>& c)
{
    std::stringstream s;
    for (const auto& e : c)
        s << e << ",";
    auto result = s.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

template <>
std::string printVectorElements(const std::vector<uint8_t>& c)
{
    std::stringstream s;
    for (const auto& e : c)
        s << static_cast<uint16_t>(e)
          << ",";  // cast to larger uint, force print as number not character
    auto result = s.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

template <>
std::string printVectorElements(const std::vector<int8_t>& c)
{
    std::stringstream s;
    for (const auto& e : c)
        s << static_cast<int16_t>(e)
          << ",";  // cast to larger int, force print as number not character
    auto result = s.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

}  // namespace pbindexdump

CppFormatter::CppFormatter(const Settings& settings) : IFormatter(settings) {}

void CppFormatter::Run()
{
    using namespace PacBio::BAM;

    const PbiRawData rawData{settings_.inputPbiFilename_};
    const PbiRawBarcodeData& barcodeData = rawData.BarcodeData();
    const PbiRawBasicData& basicData = rawData.BasicData();
    const PbiRawMappedData& mappedData = rawData.MappedData();
    const PbiRawReferenceData& referenceData = rawData.ReferenceData();

    auto version = std::string{};
    switch (rawData.Version()) {
        case PbiFile::Version_3_0_0:
            version = "PbiFile::Version_3_0_0";
            break;
        case PbiFile::Version_3_0_1:
            version = "PbiFile::Version_3_0_1";
            break;
        default:
            throw std::runtime_error("unsupported PBI version encountered");
    }

    auto fileSections = std::string{"PbiFile::BASIC"};
    if (rawData.HasBarcodeData()) fileSections += std::string{" | PbiFile::BARCODE"};
    if (rawData.HasMappedData()) fileSections += std::string{" | PbiFile::MAPPED"};
    if (rawData.HasReferenceData()) fileSections += std::string{" | PbiFile::REFERENCE"};

    std::ostringstream s;
    s << "PbiRawData rawData;" << std::endl
      << "rawData.Version(" << version << ");" << std::endl
      << "rawData.FileSections(" << fileSections << ");" << std::endl
      << "rawData.NumReads(" << rawData.NumReads() << ");" << std::endl
      << std::endl
      << "PbiRawBasicData& basicData = rawData.BasicData();" << std::endl
      << "basicData.rgId_       = {" << printVectorElements(basicData.rgId_) << "};" << std::endl
      << "basicData.qStart_     = {" << printVectorElements(basicData.qStart_) << "};" << std::endl
      << "basicData.qEnd_       = {" << printVectorElements(basicData.qEnd_) << "};" << std::endl
      << "basicData.holeNumber_ = {" << printVectorElements(basicData.holeNumber_) << "};"
      << std::endl
      << "basicData.readQual_   = {" << printVectorElements(basicData.readQual_) << "};"
      << std::endl
      << "basicData.ctxtFlag_   = {" << printVectorElements(basicData.ctxtFlag_) << "};"
      << std::endl
      << "basicData.fileOffset_ = {" << printVectorElements(basicData.fileOffset_) << "};"
      << std::endl
      << std::endl;

    if (rawData.HasBarcodeData()) {
        s << "PbiRawBarcodeData& barcodeData = rawData.BarcodeData();" << std::endl
          << "barcodeData.bcForward_ = {" << printVectorElements(barcodeData.bcForward_) << "};"
          << std::endl
          << "barcodeData.bcReverse_ = {" << printVectorElements(barcodeData.bcReverse_) << "};"
          << std::endl
          << "barcodeData.bcQual_    = {" << printVectorElements(barcodeData.bcQual_) << "};"
          << std::endl
          << std::endl;
    }

    if (rawData.HasMappedData()) {
        s << "PbiRawMappedData& mappedData = rawData.MappedData();" << std::endl
          << "mappedData.tId_       = {" << printVectorElements(mappedData.tId_) << "};"
          << std::endl
          << "mappedData.tStart_    = {" << printVectorElements(mappedData.tStart_) << "};"
          << std::endl
          << "mappedData.tEnd_      = {" << printVectorElements(mappedData.tEnd_) << "};"
          << std::endl
          << "mappedData.aStart_    = {" << printVectorElements(mappedData.aStart_) << "};"
          << std::endl
          << "mappedData.aEnd_      = {" << printVectorElements(mappedData.aEnd_) << "};"
          << std::endl
          << "mappedData.revStrand_ = {" << printVectorElements(mappedData.revStrand_) << "};"
          << std::endl
          << "mappedData.nM_        = {" << printVectorElements(mappedData.nM_) << "};" << std::endl
          << "mappedData.nMM_       = {" << printVectorElements(mappedData.nMM_) << "};"
          << std::endl
          << "mappedData.mapQV_     = {" << printVectorElements(mappedData.mapQV_) << "};"
          << std::endl
          << std::endl;
    }

    if (rawData.HasReferenceData()) {
        s << "PbiRawReferenceData& referenceData = rawData.ReferenceData();" << std::endl
          << "referenceData.entries_ = { " << std::endl
          << printCppReferenceData(referenceData) << "};" << std::endl
          << std::endl;
    }

    std::cout << s.str() << std::endl;
}

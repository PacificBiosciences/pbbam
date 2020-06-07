// Author: Derek Barnett

#include "CppFormatter.h"

#include <cstdint>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <pbbam/PbiFile.h>
#include <pbbam/PbiRawData.h>

namespace PacBio {
namespace PbIndexDump {
namespace {

std::string printReferenceData(const BAM::PbiRawReferenceData& referenceData)
{
    std::ostringstream out;
    for (const auto& entry : referenceData.entries_) {
        if (!out.str().empty()) out << ",\n";

        out << "    PbiReferenceEntry{" << entry.tId_ << "," << entry.beginRow_ << ","
            << entry.endRow_ << "}";
    }
    if (!out.str().empty()) out << '\n';
    return out.str();
}

template <typename T>
std::string printField(const std::vector<T>& c)
{
    std::ostringstream out;
    for (const auto& e : c)
        out << e << ",";
    auto result = out.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

template <>
std::string printField(const std::vector<uint8_t>& c)
{
    std::ostringstream out;
    for (const auto& e : c)
        out << static_cast<uint16_t>(e)
            << ",";  // cast to larger uint, force print as number not character
    auto result = out.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

template <>
std::string printField(const std::vector<int8_t>& c)
{
    std::ostringstream out;
    for (const auto& e : c)
        out << static_cast<int16_t>(e)
            << ",";  // cast to larger int, force print as number not character
    auto result = out.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

}  // namespace

void CppFormatter::Run(const Settings& settings)
{
    const BAM::PbiRawData rawData{settings.InputFile};
    const auto& barcodeData = rawData.BarcodeData();
    const auto& basicData = rawData.BasicData();
    const auto& mappedData = rawData.MappedData();
    const auto& referenceData = rawData.ReferenceData();

    std::string version;
    switch (rawData.Version()) {
        case BAM::PbiFile::Version_3_0_0:
            version = "PbiFile::Version_3_0_0";
            break;
        case BAM::PbiFile::Version_3_0_1:
            version = "PbiFile::Version_3_0_1";
            break;
        case BAM::PbiFile::Version_3_0_2:
            version = "PbiFile::Version_3_0_2";
            break;
        default:
            throw std::runtime_error("unsupported PBI version encountered");
    }

    std::string fileSections{"PbiFile::BASIC"};
    if (rawData.HasBarcodeData()) fileSections += std::string{" | PbiFile::BARCODE"};
    if (rawData.HasMappedData()) fileSections += std::string{" | PbiFile::MAPPED"};
    if (rawData.HasReferenceData()) fileSections += std::string{" | PbiFile::REFERENCE"};

    std::ostringstream s;
    s << "PbiRawData rawData;\n"
      << "rawData.Version(" << version << ");\n"
      << "rawData.FileSections(" << fileSections << ");\n"
      << "rawData.NumReads(" << rawData.NumReads() << ");\n"
      << '\n'
      << "PbiRawBasicData& basicData = rawData.BasicData();\n"
      << "basicData.rgId_       = {" << printField(basicData.rgId_) << "};\n"
      << "basicData.qStart_     = {" << printField(basicData.qStart_) << "};\n"
      << "basicData.qEnd_       = {" << printField(basicData.qEnd_) << "};\n"
      << "basicData.holeNumber_ = {" << printField(basicData.holeNumber_) << "};\n"
      << "basicData.readQual_   = {" << printField(basicData.readQual_) << "};\n"
      << "basicData.ctxtFlag_   = {" << printField(basicData.ctxtFlag_) << "};\n"
      << "basicData.fileOffset_ = {" << printField(basicData.fileOffset_) << "};\n";

    if (rawData.HasBarcodeData()) {
        s << '\n'
          << "PbiRawBarcodeData& barcodeData = rawData.BarcodeData();\n"
          << "barcodeData.bcForward_ = {" << printField(barcodeData.bcForward_) << "};\n"
          << "barcodeData.bcReverse_ = {" << printField(barcodeData.bcReverse_) << "};\n"
          << "barcodeData.bcQual_    = {" << printField(barcodeData.bcQual_) << "};\n";
    }

    if (rawData.HasMappedData()) {
        s << '\n'
          << "PbiRawMappedData& mappedData = rawData.MappedData();\n"
          << "mappedData.tId_       = {" << printField(mappedData.tId_) << "};\n"
          << "mappedData.tStart_    = {" << printField(mappedData.tStart_) << "};\n"
          << "mappedData.tEnd_      = {" << printField(mappedData.tEnd_) << "};\n"
          << "mappedData.aStart_    = {" << printField(mappedData.aStart_) << "};\n"
          << "mappedData.aEnd_      = {" << printField(mappedData.aEnd_) << "};\n"
          << "mappedData.revStrand_ = {" << printField(mappedData.revStrand_) << "};\n"
          << "mappedData.nM_        = {" << printField(mappedData.nM_) << "};\n"
          << "mappedData.nMM_       = {" << printField(mappedData.nMM_) << "};\n"
          << "mappedData.mapQV_     = {" << printField(mappedData.mapQV_) << "};\n";
    }

    if (rawData.HasReferenceData()) {
        s << '\n'
          << "PbiRawReferenceData& referenceData = rawData.ReferenceData();\n"
          << "referenceData.entries_ = { \n"
          << printReferenceData(referenceData) << "};\n";
    }

    std::cout << s.str();
}

}  // namespace PbIndexDump
}  // namespace PacBio

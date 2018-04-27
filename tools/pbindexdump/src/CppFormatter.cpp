// Author: Derek Barnett

#include "CppFormatter.h"
#include <pbbam/PbiRawData.h>

#include <cstdint>
#include <iostream>
#include <sstream>

using namespace pbindexdump;
using namespace std;

namespace pbindexdump {

static string printCppReferenceData(const PacBio::BAM::PbiRawReferenceData& referenceData)
{
    auto result = string{""};
    for (const PacBio::BAM::PbiReferenceEntry& entry : referenceData.entries_) {
        if (!result.empty()) result.append(",\n");
        result.append(string{"    PbiReferenceEntry{"} + to_string(entry.tId_) + "," +
                      to_string(entry.beginRow_) + "," + to_string(entry.endRow_) + string{"}"});
    }
    if (!result.empty()) result.append("\n");
    return result;
}

template <typename T>
string printVectorElements(const std::vector<T>& c)
{
    stringstream s;
    for (const auto& e : c)
        s << e << ",";
    auto result = s.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

template <>
string printVectorElements(const std::vector<uint8_t>& c)
{
    stringstream s;
    for (const auto& e : c)
        s << static_cast<uint16_t>(e)
          << ",";  // cast to larger uint, force print as number not character
    auto result = s.str();
    if (!result.empty()) result.pop_back();  // remove final comma
    return result;
}

template <>
string printVectorElements(const std::vector<int8_t>& c)
{
    stringstream s;
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

    auto version = string{};
    switch (rawData.Version()) {
        case PbiFile::Version_3_0_0:
            version = "PbiFile::Version_3_0_0";
            break;
        case PbiFile::Version_3_0_1:
            version = "PbiFile::Version_3_0_1";
            break;
        default:
            throw runtime_error("unsupported PBI version encountered");
    }

    auto fileSections = string{"PbiFile::BASIC"};
    if (rawData.HasBarcodeData()) fileSections += string{" | PbiFile::BARCODE"};
    if (rawData.HasMappedData()) fileSections += string{" | PbiFile::MAPPED"};
    if (rawData.HasReferenceData()) fileSections += string{" | PbiFile::REFERENCE"};

    stringstream s;
    s << "PbiRawData rawData;" << endl
      << "rawData.Version(" << version << ");" << endl
      << "rawData.FileSections(" << fileSections << ");" << endl
      << "rawData.NumReads(" << rawData.NumReads() << ");" << endl
      << endl
      << "PbiRawBasicData& basicData = rawData.BasicData();" << endl
      << "basicData.rgId_       = {" << printVectorElements(basicData.rgId_) << "};" << endl
      << "basicData.qStart_     = {" << printVectorElements(basicData.qStart_) << "};" << endl
      << "basicData.qEnd_       = {" << printVectorElements(basicData.qEnd_) << "};" << endl
      << "basicData.holeNumber_ = {" << printVectorElements(basicData.holeNumber_) << "};" << endl
      << "basicData.readQual_   = {" << printVectorElements(basicData.readQual_) << "};" << endl
      << "basicData.ctxtFlag_   = {" << printVectorElements(basicData.ctxtFlag_) << "};" << endl
      << "basicData.fileOffset_ = {" << printVectorElements(basicData.fileOffset_) << "};" << endl
      << endl;

    if (rawData.HasBarcodeData()) {
        s << "PbiRawBarcodeData& barcodeData = rawData.BarcodeData();" << endl
          << "barcodeData.bcForward_ = {" << printVectorElements(barcodeData.bcForward_) << "};"
          << endl
          << "barcodeData.bcReverse_ = {" << printVectorElements(barcodeData.bcReverse_) << "};"
          << endl
          << "barcodeData.bcQual_    = {" << printVectorElements(barcodeData.bcQual_) << "};"
          << endl
          << endl;
    }

    if (rawData.HasMappedData()) {
        s << "PbiRawMappedData& mappedData = rawData.MappedData();" << endl
          << "mappedData.tId_       = {" << printVectorElements(mappedData.tId_) << "};" << endl
          << "mappedData.tStart_    = {" << printVectorElements(mappedData.tStart_) << "};" << endl
          << "mappedData.tEnd_      = {" << printVectorElements(mappedData.tEnd_) << "};" << endl
          << "mappedData.aStart_    = {" << printVectorElements(mappedData.aStart_) << "};" << endl
          << "mappedData.aEnd_      = {" << printVectorElements(mappedData.aEnd_) << "};" << endl
          << "mappedData.revStrand_ = {" << printVectorElements(mappedData.revStrand_) << "};"
          << endl
          << "mappedData.nM_        = {" << printVectorElements(mappedData.nM_) << "};" << endl
          << "mappedData.nMM_       = {" << printVectorElements(mappedData.nMM_) << "};" << endl
          << "mappedData.mapQV_     = {" << printVectorElements(mappedData.mapQV_) << "};" << endl
          << endl;
    }

    if (rawData.HasReferenceData()) {
        s << "PbiRawReferenceData& referenceData = rawData.ReferenceData();" << endl
          << "referenceData.entries_ = { " << endl
          << printCppReferenceData(referenceData) << "};" << endl
          << endl;
    }

    cout << s.str() << endl;
}

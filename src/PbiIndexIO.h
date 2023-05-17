#ifndef PBBAM_PBIINDEXIO_H
#define PBBAM_PBIINDEXIO_H

#include <pbbam/Config.h>

#include <pbbam/BamFile.h>
#include <pbbam/DataSet.h>
#include <pbbam/PbiFile.h>
#include <pbbam/PbiRawData.h>

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

class PbiIndexIO
{
public:
    // top-level entry points
    static void LoadFromFile(PbiRawData& rawData, const std::string& filename);
    static void LoadFromDataSet(PbiRawData& aggregateData, const DataSet& dataset);
    static void Save(const PbiRawData& rawData, const std::string& filename);

    // per-component load
    static void LoadBarcodeData(PbiRawBarcodeData& barcodeData, uint32_t numReads, BGZF* fp);
    static void LoadHeader(PbiRawData& index, BGZF* fp);
    static void LoadMappedData(PbiRawMappedData& mappedData, uint32_t numReads, BGZF* fp);
    static void LoadReferenceData(PbiRawReferenceData& referenceData, BGZF* fp);
    static void LoadBasicData(PbiRawBasicData& basicData, uint32_t numReads, BGZF* fp);

    // per-data-field load
    template <typename T>
    static void LoadBgzfVector(BGZF* fp, std::vector<T>& data, uint32_t numReads);

    // per-component write
    static void WriteBarcodeData(const PbiRawBarcodeData& barcodeData, uint32_t numReads, BGZF* fp);
    static void WriteHeader(const PbiRawData& index, BGZF* fp);
    static void WriteMappedData(const PbiRawMappedData& mappedData, uint32_t numReads, BGZF* fp);
    static void WriteReferenceData(const PbiRawReferenceData& referenceData, BGZF* fp);
    static void WriteBasicData(const PbiRawBasicData& subreadData, uint32_t numReads, BGZF* fp);

    // per-data-field write
    template <typename T>
    static void WriteBgzfVector(BGZF* fp, const std::vector<T>& data);

private:
    // helper functions
    template <typename T>
    static void SwapEndianness(std::vector<T>& data);
};

template <typename T>
void PbiIndexIO::LoadBgzfVector(BGZF* fp, std::vector<T>& data, const uint32_t numReads)
{
    assert(fp);
    data.resize(numReads);
    auto ret = bgzf_read(fp, &data[0], numReads * sizeof(T));
    if (fp->is_be) {
        SwapEndianness(data);
    }
    std::ignore = ret;
}

template <typename T>
void PbiIndexIO::SwapEndianness(std::vector<T>& data)
{
    const auto elementSize = sizeof(T);
    const auto numReads = data.size();
    switch (elementSize) {
        case 1:
            break;  // no swapping necessary
        case 2:
            for (std::size_t i = 0; i < numReads; ++i) {
                ed_swap_2p(&data[i]);
            }
            break;
        case 4:
            for (std::size_t i = 0; i < numReads; ++i) {
                ed_swap_4p(&data[i]);
            }
            break;
        case 8:
            for (std::size_t i = 0; i < numReads; ++i) {
                ed_swap_8p(&data[i]);
            }
            break;
        default:
            throw std::runtime_error{"unsupported element size"};
    }
}

template <typename T>
void PbiIndexIO::WriteBgzfVector(BGZF* fp, const std::vector<T>& data)
{
    assert(fp);
    std::vector<T> output = data;
    if (fp->is_be) {
        SwapEndianness(output);
    }
    auto ret = bgzf_write(fp, &output[0], data.size() * sizeof(T));
    std::ignore = ret;
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_PBIINDEXIO_H

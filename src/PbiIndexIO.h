// Author: Derek Barnett

#ifndef PBIINDEXIO_H
#define PBIINDEXIO_H

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include "MemoryUtils.h"
#include "pbbam/BamFile.h"
#include "pbbam/DataSet.h"
#include "pbbam/PbiFile.h"
#include "pbbam/PbiRawData.h"
#include "pbbam/Unused.h"

namespace PacBio {
namespace BAM {

struct PbiHeader
{
    uint32_t numReads = 0;
    PbiFile::VersionEnum version = PbiFile::CurrentVersion;
    PbiFile::Sections sections = PbiFile::ALL;
    int64_t firstRecordOffset;  // PBI, not BAM. only used internally by PbiIndexIO
};

class PbiIndexIO
{
public:
    // special dataset-handler
    static PbiRawData LoadFromDataSet(const DataSet& dataset);

public:
    explicit PbiIndexIO(const std::string& pbiFilename);
    PbiIndexIO(const std::string& pbiFilename, const std::set<PbiFile::Field>& fields);

    PbiRawData Load();
    const PbiHeader& Header() const;

private:
    std::string pbiFilename_;  // empty if from dataset
    std::set<PbiFile::Field> fields_;
    std::unique_ptr<BGZF, HtslibBgzfDeleter> fp_;
    PbiHeader header_;

    // dummy buffer for skipped fields
    std::vector<int64_t> temp_;

    void LoadBarcodeData(PbiRawData& data);
    void LoadBasicData(PbiRawData& data);
    void LoadHeader();
    void LoadMappedData(PbiRawData& data);
    void LoadReferenceData(PbiRawData& data);
    void Open(const std::string& filename);

    template <typename T>
    void MaybeSaveField(std::vector<T>& dst, const PbiFile::Field field);

    template <typename T>
    void SaveField(std::vector<T>& dst);

    template <size_t ElementSize>
    void SkipField();

    static void AggregateDataSet(PbiRawData& aggregateData, const DataSet& dataset);

    // // per-data-field load
    template <typename T>
    static void LoadBgzfVector(BGZF* fp, std::vector<T>& data, const uint32_t numReads,
                               const bool maybeSwapEndian = true);

private:
    // helper functions
    template <typename T>
    static void SwapEndianness(std::vector<T>& data);
};

template <typename T>
inline void PbiIndexIO::LoadBgzfVector(BGZF* fp, std::vector<T>& data, const uint32_t numReads,
                                       const bool maybeSwapEndian)
{
    assert(fp);
    data.resize(numReads);
    auto ret = bgzf_read(fp, &data[0], numReads * sizeof(T));
    if (fp->is_be && maybeSwapEndian) SwapEndianness(data);
    UNUSED(ret);
}

template <typename T>
inline void PbiIndexIO::SwapEndianness(std::vector<T>& data)
{
    const auto elementSize = sizeof(T);
    const auto numReads = data.size();
    switch (elementSize) {
        case 1:
            break;  // no swapping necessary
        case 2:
            for (size_t i = 0; i < numReads; ++i)
                ed_swap_2p(&data[i]);
            break;
        case 4:
            for (size_t i = 0; i < numReads; ++i)
                ed_swap_4p(&data[i]);
            break;
        case 8:
            for (size_t i = 0; i < numReads; ++i)
                ed_swap_8p(&data[i]);
            break;
        default:
            throw std::runtime_error{"unsupported element size"};
    }
}

}  // namespace BAM
}  // namespace PacBio

#endif  // PBIINDEXIO_H

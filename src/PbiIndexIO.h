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
//
// Author: Derek Barnett

#ifndef PBIINDEXIO_H
#define PBIINDEXIO_H

#include "pbbam/BamFile.h"
#include "pbbam/DataSet.h"
#include "pbbam/PbiFile.h"
#include "pbbam/PbiRawData.h"
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <memory>
#include <string>
#include <vector>
#include <cassert>

namespace PacBio {
namespace BAM {
namespace internal {

class PbiIndexIO
{
public:
    // top-level entry points
    static PbiRawData Load(const std::string& filename);
    static void Load(PbiRawData& rawData, const std::string& filename);
    static void LoadFromDataSet(PbiRawData& aggregateData, const DataSet& dataset);
    static void Save(const PbiRawData& rawData, const std::string& filename);

public:
    // per-component load
    static void LoadBarcodeData(PbiRawBarcodeData& barcodeData,
                                const uint32_t numReads,
                                BGZF* fp);
    static void LoadHeader(PbiRawData& index,
                           BGZF* fp);
    static void LoadMappedData(PbiRawMappedData& mappedData,
                               const uint32_t numReads,
                               BGZF* fp);
    static void LoadReferenceData(PbiRawReferenceData& referenceData,
                                  BGZF* fp);
    static void LoadBasicData(PbiRawBasicData& basicData,
                              const uint32_t numReads,
                              BGZF* fp);

    // per-data-field load
    template<typename T>
    static void LoadBgzfVector(BGZF* fp,
                               std::vector<T>& data,
                               const uint32_t numReads);

public:
    // per-component write
    static void WriteBarcodeData(const PbiRawBarcodeData& barcodeData,
                                 const uint32_t numReads,
                                 BGZF* fp);
    static void WriteHeader(const PbiRawData& index,
                            BGZF* fp);
    static void WriteMappedData(const PbiRawMappedData& mappedData,
                                const uint32_t numReads,
                                BGZF* fp);
    static void WriteReferenceData(const PbiRawReferenceData& referenceData,
                                   BGZF* fp);
    static void WriteBasicData(const PbiRawBasicData& subreadData,
                                 const uint32_t numReads,
                                 BGZF* fp);

    // per-data-field write
    template<typename T>
    static void WriteBgzfVector(BGZF* fp,
                                const std::vector<T>& data);

private:
    // helper functions
    template<typename T>
    static void SwapEndianness(std::vector<T>& data);
};

template<typename T>
inline void PbiIndexIO::LoadBgzfVector(BGZF* fp,
                                       std::vector<T>& data,
                                       const uint32_t numReads)
{
    assert(fp);
    data.resize(numReads);
    bgzf_read(fp, &data[0], numReads*sizeof(T));
    if (fp->is_be)
        SwapEndianness(data);
}

template<typename T>
inline void PbiIndexIO::SwapEndianness(std::vector<T>& data)
{
    const size_t elementSize = sizeof(T);
    const size_t numReads = data.size();
    switch (elementSize) {
        case 1 : break; // no swapping necessary
        case 2 :
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
            throw std::runtime_error("unsupported element size");
    }
}

template<typename T>
inline void PbiIndexIO::WriteBgzfVector(BGZF* fp,
                                        const std::vector<T>& data)
{
    assert(fp);
    std::vector<T> output = data;
    if (fp->is_be)
        SwapEndianness(output);
    bgzf_write(fp, &output[0], data.size()*sizeof(T));
}

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // PBIINDEXIO_H

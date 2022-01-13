#include "PbbamInternalConfig.h"

#include "PbiIndexIO.h"

#include <pbbam/BamRecord.h>
#include <pbbam/Deleters.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiBuilder.h>
#include "ErrnoReason.h"

#include <pbcopper/utility/MoveAppend.h>

#include <boost/algorithm/string.hpp>

#include <array>
#include <optional>
#include <sstream>
#include <stdexcept>

#include <cassert>

namespace PacBio {
namespace BAM {
namespace {

void CheckContainer(const std::string& container, const size_t expected, const size_t observed)
{
    if (observed != expected) {
        std::ostringstream msg;
        msg << "[pbbam] PBI index I/O ERROR: expected " << expected << " records in " << container
            << " field, but found " << observed << " instead";
        throw std::runtime_error{msg.str()};
    }
}

void CheckExpectedSize(const PbiRawBarcodeData& barcodeData, const size_t numReads)
{
    CheckContainer("BarcodeData.bc_forward", numReads, barcodeData.bcForward_.size());
    CheckContainer("BarcodeData.bc_reverse", numReads, barcodeData.bcReverse_.size());
    CheckContainer("BarcodeData.bc_qual", numReads, barcodeData.bcReverse_.size());
}

void CheckExpectedSize(const PbiRawBasicData& basicData, const size_t numReads)
{
    CheckContainer("BasicData.rgId", numReads, basicData.rgId_.size());
    CheckContainer("BasicData.qStart", numReads, basicData.qStart_.size());
    CheckContainer("BasicData.qEnd", numReads, basicData.qEnd_.size());
    CheckContainer("BasicData.holeNumber", numReads, basicData.holeNumber_.size());
    CheckContainer("BasicData.readQual", numReads, basicData.readQual_.size());
    CheckContainer("BasicData.ctxt_flag", numReads, basicData.ctxtFlag_.size());
    CheckContainer("BasicData.fileOffset", numReads, basicData.fileOffset_.size());
}

void CheckExpectedSize(const PbiRawMappedData& mappedData, const size_t numReads)
{
    CheckContainer("MappedData.tId", numReads, mappedData.tId_.size());
    CheckContainer("MappedData.tStart", numReads, mappedData.tStart_.size());
    CheckContainer("MappedData.tEnd", numReads, mappedData.tEnd_.size());
    CheckContainer("MappedData.aStart", numReads, mappedData.aStart_.size());
    CheckContainer("MappedData.aEnd", numReads, mappedData.aEnd_.size());
    CheckContainer("MappedData.revStrand", numReads, mappedData.revStrand_.size());
    CheckContainer("MappedData.nM", numReads, mappedData.nM_.size());
    CheckContainer("MappedData.nMM", numReads, mappedData.nMM_.size());
    CheckContainer("MappedData.mapQV", numReads, mappedData.mapQV_.size());

    if (mappedData.hasIndelOps_) {
        CheckContainer("MappedData.nInsOps", numReads, mappedData.nInsOps_.size());
        CheckContainer("MappedData.nDelOps", numReads, mappedData.nDelOps_.size());
    }
}

}  // namespace

void PbiIndexIO::LoadFromFile(PbiRawData& rawData, const std::string& filename)
{
    // open file for reading
    if (!boost::algorithm::iends_with(filename, ".pbi")) {
        std::ostringstream msg;
        msg << "[pbbam] PBI index I/O ERROR: unsupported file extension:\n"
            << "  file: " << filename;
        throw std::runtime_error{msg.str()};
    }

    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf(bgzf_open(filename.c_str(), "rb"));
    auto* fp = bgzf.get();
    if (fp == nullptr) {
        std::ostringstream msg;
        msg << "[pbbam] PBI index I/O ERROR: could not open file for reading:\n"
            << "  file: " << filename;
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }

    // load data
    LoadHeader(rawData, fp);
    const auto numReads = rawData.NumReads();
    if (numReads > 0) {
        LoadBasicData(rawData.BasicData(), numReads, fp);
        if (rawData.HasMappedData()) {
            LoadMappedData(rawData.MappedData(), numReads, fp);
        }
        if (rawData.HasReferenceData()) {
            LoadReferenceData(rawData.ReferenceData(), fp);
        }
        if (rawData.HasBarcodeData()) {
            LoadBarcodeData(rawData.BarcodeData(), numReads, fp);
        }
    }
}

void PbiIndexIO::LoadFromDataSet(PbiRawData& aggregateData, const DataSet& dataset)
{
    aggregateData.NumReads(0);
    aggregateData.FileSections(PbiFile::BASIC | PbiFile::MAPPED | PbiFile::BARCODE);

    // Some GCC configurations give false-positive warnings against using uninitialized
    // std::optional here, hence the 'old-fashioned' bool flag.
    bool isSet = false;
    PbiFile::VersionEnum aggregateVersion = PbiFile::CurrentVersion;
    const auto compatibleVersion = [&](PbiFile::VersionEnum next) {
        if (!isSet) {
            aggregateVersion = next;
            isSet = true;
            return true;
        } else {
            return (aggregateVersion < PbiFile::Version_4_0_0 && next < PbiFile::Version_4_0_0) ||
                   (aggregateVersion >= PbiFile::Version_4_0_0 && next >= PbiFile::Version_4_0_0);
        }
    };

    const auto bamFiles = dataset.BamFiles();
    uint16_t fileNumber = 0;
    for (const auto& bamFile : bamFiles) {
        PbiRawData currentPbi{bamFile.PacBioIndexFilename()};
        if (!compatibleVersion(currentPbi.Version())) {
            throw std::runtime_error{
                "[pbbam] PBI index I/O ERROR: dataset contains incompatible PBI index versions. "
                "Please rerun BAM files through 'pbindex' to ensure compatibility."};
        }

        // read count
        const auto currentPbiCount = currentPbi.NumReads();
        aggregateData.NumReads(aggregateData.NumReads() + currentPbiCount);

        // BasicData
        auto& aggregateBasicData = aggregateData.BasicData();
        auto& currentBasicData = currentPbi.BasicData();
        Utility::MoveAppend(std::move(currentBasicData.rgId_), aggregateBasicData.rgId_);
        Utility::MoveAppend(std::move(currentBasicData.qStart_), aggregateBasicData.qStart_);
        Utility::MoveAppend(std::move(currentBasicData.qEnd_), aggregateBasicData.qEnd_);
        Utility::MoveAppend(std::move(currentBasicData.holeNumber_),
                            aggregateBasicData.holeNumber_);
        Utility::MoveAppend(std::move(currentBasicData.readQual_), aggregateBasicData.readQual_);
        Utility::MoveAppend(std::move(currentBasicData.ctxtFlag_), aggregateBasicData.ctxtFlag_);
        Utility::MoveAppend(std::move(currentBasicData.fileOffset_),
                            aggregateBasicData.fileOffset_);
        Utility::MoveAppend(std::vector<uint16_t>(currentPbiCount, fileNumber),
                            aggregateBasicData.fileNumber_);

        // BarcodeData
        auto& aggregateBarcodeData = aggregateData.BarcodeData();
        if (currentPbi.HasBarcodeData()) {
            auto& currentBarcodeData = currentPbi.BarcodeData();
            Utility::MoveAppend(std::move(currentBarcodeData.bcForward_),
                                aggregateBarcodeData.bcForward_);
            Utility::MoveAppend(std::move(currentBarcodeData.bcReverse_),
                                aggregateBarcodeData.bcReverse_);
            Utility::MoveAppend(std::move(currentBarcodeData.bcQual_),
                                aggregateBarcodeData.bcQual_);
        } else {
            Utility::MoveAppend(std::vector<int16_t>(currentPbiCount, -1),
                                aggregateBarcodeData.bcForward_);
            Utility::MoveAppend(std::vector<int16_t>(currentPbiCount, -1),
                                aggregateBarcodeData.bcReverse_);
            Utility::MoveAppend(std::vector<int8_t>(currentPbiCount, -1),
                                aggregateBarcodeData.bcQual_);
        }

        // MappedData
        auto& aggregateMappedData = aggregateData.MappedData();
        if (currentPbi.HasMappedData()) {
            auto& currentMappedData = currentPbi.MappedData();
            Utility::MoveAppend(std::move(currentMappedData.tId_), aggregateMappedData.tId_);
            Utility::MoveAppend(std::move(currentMappedData.tStart_), aggregateMappedData.tStart_);
            Utility::MoveAppend(std::move(currentMappedData.tEnd_), aggregateMappedData.tEnd_);
            Utility::MoveAppend(std::move(currentMappedData.aStart_), aggregateMappedData.aStart_);
            Utility::MoveAppend(std::move(currentMappedData.aEnd_), aggregateMappedData.aEnd_);
            Utility::MoveAppend(std::move(currentMappedData.revStrand_),
                                aggregateMappedData.revStrand_);
            Utility::MoveAppend(std::move(currentMappedData.nM_), aggregateMappedData.nM_);
            Utility::MoveAppend(std::move(currentMappedData.nMM_), aggregateMappedData.nMM_);
            Utility::MoveAppend(std::move(currentMappedData.mapQV_), aggregateMappedData.mapQV_);
            if (aggregateVersion >= PbiFile::Version_4_0_0) {
                Utility::MoveAppend(std::move(currentMappedData.nInsOps_),
                                    aggregateMappedData.nInsOps_);
                Utility::MoveAppend(std::move(currentMappedData.nDelOps_),
                                    aggregateMappedData.nDelOps_);
            }

        } else {
            Utility::MoveAppend(std::vector<int32_t>(currentPbiCount, -1),
                                aggregateMappedData.tId_);
            Utility::MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                                aggregateMappedData.tStart_);
            Utility::MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                                aggregateMappedData.tEnd_);
            Utility::MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                                aggregateMappedData.aStart_);
            Utility::MoveAppend(std::vector<uint32_t>(currentPbiCount, UnmappedPosition),
                                aggregateMappedData.aEnd_);
            Utility::MoveAppend(std::vector<uint8_t>(currentPbiCount, 0),
                                aggregateMappedData.revStrand_);
            Utility::MoveAppend(std::vector<uint32_t>(currentPbiCount, 0), aggregateMappedData.nM_);
            Utility::MoveAppend(std::vector<uint32_t>(currentPbiCount, 0),
                                aggregateMappedData.nMM_);
            Utility::MoveAppend(std::vector<uint8_t>(currentPbiCount, 255),
                                aggregateMappedData.mapQV_);
            if (aggregateVersion >= PbiFile::Version_4_0_0) {
                Utility::MoveAppend(std::vector<uint8_t>(currentPbiCount, 0),
                                    aggregateMappedData.nInsOps_);
                Utility::MoveAppend(std::vector<uint8_t>(currentPbiCount, 0),
                                    aggregateMappedData.nDelOps_);
            }
        }

        ++fileNumber;
    }

    aggregateData.Version(aggregateVersion);
}

void PbiIndexIO::LoadBarcodeData(PbiRawBarcodeData& barcodeData, const uint32_t numReads, BGZF* fp)
{
    // read from file
    LoadBgzfVector(fp, barcodeData.bcForward_, numReads);
    LoadBgzfVector(fp, barcodeData.bcReverse_, numReads);
    LoadBgzfVector(fp, barcodeData.bcQual_, numReads);

    // validate
    CheckExpectedSize(barcodeData, numReads);
}

void PbiIndexIO::LoadHeader(PbiRawData& index, BGZF* fp)
{
    // 'magic' string
    char magic[4];
    auto bytesRead = bgzf_read(fp, magic, 4);
    if (bytesRead != 4 || strncmp(magic, "PBI\1", 4)) {
        throw std::runtime_error{
            "[pbbam] PBI index I/O ERROR: expected PBI file, found unknown format instead"};
    }

    // version, pbi_flags, & n_reads
    uint32_t version;
    uint16_t sections;
    uint32_t numReads;
    bytesRead = bgzf_read(fp, &version, sizeof(version));
    bytesRead = bgzf_read(fp, &sections, sizeof(sections));
    bytesRead = bgzf_read(fp, &numReads, sizeof(numReads));
    if (fp->is_be) {
        version = ed_swap_4(version);
        sections = ed_swap_2(sections);
        numReads = ed_swap_4(numReads);
    }

    index.Version(static_cast<PbiFile::VersionEnum>(version));
    index.FileSections(sections);
    index.NumReads(numReads);

    if (static_cast<PbiFile::VersionEnum>(version) < PbiFile::Version_4_0_0) {
        index.MappedData().hasIndelOps_ = false;
    }

    // skip reserved section
    size_t reservedLength = 18;
    // adjust depending on version
    char reserved[18];
    bytesRead = bgzf_read(fp, &reserved, reservedLength);
}

void PbiIndexIO::LoadMappedData(PbiRawMappedData& mappedData, const uint32_t numReads, BGZF* fp)
{
    // read from file
    LoadBgzfVector(fp, mappedData.tId_, numReads);
    LoadBgzfVector(fp, mappedData.tStart_, numReads);
    LoadBgzfVector(fp, mappedData.tEnd_, numReads);
    LoadBgzfVector(fp, mappedData.aStart_, numReads);
    LoadBgzfVector(fp, mappedData.aEnd_, numReads);
    LoadBgzfVector(fp, mappedData.revStrand_, numReads);
    LoadBgzfVector(fp, mappedData.nM_, numReads);
    LoadBgzfVector(fp, mappedData.nMM_, numReads);
    LoadBgzfVector(fp, mappedData.mapQV_, numReads);

    if (mappedData.hasIndelOps_) {
        LoadBgzfVector(fp, mappedData.nInsOps_, numReads);
        LoadBgzfVector(fp, mappedData.nDelOps_, numReads);
    }

    // validate
    CheckExpectedSize(mappedData, numReads);
}

void PbiIndexIO::LoadReferenceData(PbiRawReferenceData& referenceData, BGZF* fp)
{
    assert(sizeof(PbiReferenceEntry::ID) == 4);
    assert(sizeof(PbiReferenceEntry::Row) == 4);

    // num refs
    uint32_t numRefs;
    auto ret = bgzf_read(fp, &numRefs, 4);
    if (fp->is_be) {
        numRefs = ed_swap_4(numRefs);
    }

    // reference entries
    referenceData.entries_.clear();
    referenceData.entries_.resize(numRefs);
    for (auto& entry : referenceData.entries_) {
        ret = bgzf_read(fp, &entry.tId_, 4);
        ret = bgzf_read(fp, &entry.beginRow_, 4);
        ret = bgzf_read(fp, &entry.endRow_, 4);
        if (fp->is_be) {
            entry.tId_ = ed_swap_4(entry.tId_);
            entry.beginRow_ = ed_swap_4(entry.beginRow_);
            entry.endRow_ = ed_swap_4(entry.endRow_);
        }
    }
    std::ignore = ret;
}

void PbiIndexIO::LoadBasicData(PbiRawBasicData& basicData, const uint32_t numReads, BGZF* fp)
{
    // read from file
    LoadBgzfVector(fp, basicData.rgId_, numReads);
    LoadBgzfVector(fp, basicData.qStart_, numReads);
    LoadBgzfVector(fp, basicData.qEnd_, numReads);
    LoadBgzfVector(fp, basicData.holeNumber_, numReads);
    LoadBgzfVector(fp, basicData.readQual_, numReads);
    LoadBgzfVector(fp, basicData.ctxtFlag_, numReads);
    LoadBgzfVector(fp, basicData.fileOffset_, numReads);

    // validate
    CheckExpectedSize(basicData, numReads);
}

void PbiIndexIO::Save(const PbiRawData& index, const std::string& filename)
{
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf(bgzf_open(filename.c_str(), "wb"));
    auto* fp = bgzf.get();
    if (fp == nullptr) {
        std::ostringstream msg;
        msg << "[pbbam] PBI index I/O ERROR: could not open file for writing:\n"
            << "  file: " << filename;
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }

    WriteHeader(index, fp);
    const auto numReads = index.NumReads();
    if (numReads > 0) {
        WriteBasicData(index.BasicData(), numReads, fp);

        if (index.HasMappedData()) {
            WriteMappedData(index.MappedData(), numReads, fp);
        }
        if (index.HasReferenceData()) {
            WriteReferenceData(index.ReferenceData(), fp);
        }
        if (index.HasBarcodeData()) {
            WriteBarcodeData(index.BarcodeData(), numReads, fp);
        }
    }
}

void PbiIndexIO::WriteBarcodeData(const PbiRawBarcodeData& barcodeData, const uint32_t numReads,
                                  BGZF* fp)
{
    // validate
    CheckExpectedSize(barcodeData, numReads);

    // write to file
    WriteBgzfVector(fp, barcodeData.bcForward_);
    WriteBgzfVector(fp, barcodeData.bcReverse_);
    WriteBgzfVector(fp, barcodeData.bcQual_);
}

void PbiIndexIO::WriteHeader(const PbiRawData& index, BGZF* fp)
{
    // 'magic' string
    constexpr static const std::array<char, 4> magic{{'P', 'B', 'I', '\1'}};
    auto ret = bgzf_write(fp, magic.data(), 4);

    // version, pbi_flags, & n_reads
    auto version = static_cast<uint32_t>(index.Version());
    uint16_t pbi_flags = index.FileSections();
    auto numReads = static_cast<uint16_t>(index.NumReads());
    if (fp->is_be) {
        version = ed_swap_4(version);
        pbi_flags = ed_swap_2(pbi_flags);
        numReads = ed_swap_4(numReads);
    }
    ret = bgzf_write(fp, &version, 4);
    ret = bgzf_write(fp, &pbi_flags, 2);
    ret = bgzf_write(fp, &numReads, 4);

    // reserved space
    char reserved[18];
    memset(reserved, 0, 18);
    ret = bgzf_write(fp, reserved, 18);
    std::ignore = ret;
}

void PbiIndexIO::WriteMappedData(const PbiRawMappedData& mappedData, const uint32_t numReads,
                                 BGZF* fp)
{
    // validate
    CheckExpectedSize(mappedData, numReads);

    // write to file
    WriteBgzfVector(fp, mappedData.tId_);
    WriteBgzfVector(fp, mappedData.tStart_);
    WriteBgzfVector(fp, mappedData.tEnd_);
    WriteBgzfVector(fp, mappedData.aStart_);
    WriteBgzfVector(fp, mappedData.aEnd_);
    WriteBgzfVector(fp, mappedData.revStrand_);
    WriteBgzfVector(fp, mappedData.nM_);
    WriteBgzfVector(fp, mappedData.nMM_);
    WriteBgzfVector(fp, mappedData.mapQV_);

    if (mappedData.hasIndelOps_) {
        WriteBgzfVector(fp, mappedData.nInsOps_);
        WriteBgzfVector(fp, mappedData.nDelOps_);
    }
}

void PbiIndexIO::WriteReferenceData(const PbiRawReferenceData& referenceData, BGZF* fp)
{
    // num_refs
    auto numRefs = referenceData.entries_.size();
    if (fp->is_be) {
        numRefs = ed_swap_4(numRefs);
    }
    auto ret = bgzf_write(fp, &numRefs, 4);

    // reference entries
    for (const auto& entry : referenceData.entries_) {
        auto tId = entry.tId_;
        auto beginRow = entry.beginRow_;
        auto endRow = entry.endRow_;
        if (fp->is_be) {
            tId = ed_swap_4(tId);
            beginRow = ed_swap_4(beginRow);
            endRow = ed_swap_4(endRow);
        }
        ret = bgzf_write(fp, &tId, 4);
        ret = bgzf_write(fp, &beginRow, 4);
        ret = bgzf_write(fp, &endRow, 4);
    }
    std::ignore = ret;
}

void PbiIndexIO::WriteBasicData(const PbiRawBasicData& basicData, const uint32_t numReads, BGZF* fp)
{
    // validate
    CheckExpectedSize(basicData, numReads);

    // write to file
    WriteBgzfVector(fp, basicData.rgId_);
    WriteBgzfVector(fp, basicData.qStart_);
    WriteBgzfVector(fp, basicData.qEnd_);
    WriteBgzfVector(fp, basicData.holeNumber_);
    WriteBgzfVector(fp, basicData.readQual_);
    WriteBgzfVector(fp, basicData.ctxtFlag_);
    WriteBgzfVector(fp, basicData.fileOffset_);
}

}  // namespace BAM
}  // namespace PacBio

// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamWriter.h"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <thread>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "Autovalidate.h"
#include "FileProducer.h"
#include "MemoryUtils.h"
#include "pbbam/BamFile.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/Unused.h"
#include "pbbam/Validator.h"

namespace PacBio {
namespace BAM {

class BamWriter::BamWriterPrivate
{
public:
    BamWriterPrivate(const std::string& filename, const std::shared_ptr<bam_hdr_t> rawHeader,
                     const BamWriter::CompressionLevel compressionLevel, const size_t numThreads,
                     const BamWriter::BinCalculationMode binCalculationMode, const bool useTempFile)
        : calculateBins_{binCalculationMode == BamWriter::BinCalculation_ON}, header_{rawHeader}
    {
        if (!header_) throw std::runtime_error{"BamWriter: null header provided for: " + filename};

        if (useTempFile) fileProducer_ = std::make_unique<FileProducer>(filename);

        // open file
        const auto usingFilename = (fileProducer_ ? fileProducer_->TempFilename() : filename);
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel));
        file_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!file_)
            throw std::runtime_error{"BamWriter: could not open BAM file for writing: " +
                                     usingFilename};

        // if no explicit thread count given, attempt built-in check
        size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) actualNumThreads = 1;
        }

        // if multithreading requested, enable it
        if (actualNumThreads > 1) hts_set_threads(file_.get(), actualNumThreads);

        // write header
        const auto ret = sam_hdr_write(file_.get(), header_.get());
        if (ret != 0)
            throw std::runtime_error{"BamWriter: could not write header for file: " +
                                     usingFilename};
    }

    void Write(const BamRecord& record)
    {
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif

        const auto rawRecord = BamRecordMemory::GetRawData(record);

        // (probably) store bins
        // min_shift=14 & n_lvls=5 are BAM "magic numbers"
        if (calculateBins_)
            rawRecord->core.bin =
                hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file
        const auto ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw std::runtime_error{"BamWriter: could not write record to file"};
    }

    void Write(const BamRecord& record, int64_t* vOffset)
    {
        BGZF* bgzf = file_.get()->fp.bgzf;
        assert(bgzf);
        assert(vOffset);

        // ensure offsets up-to-date
        const auto ret = bgzf_flush(bgzf);
        UNUSED(ret);

        // capture virtual offset where we’re about to write
        const auto rawTell = htell(bgzf->fp);
        const auto length = bgzf->block_offset;
        *vOffset = (rawTell << 16) | length;

        // now write data
        Write(record);
    }

    void Write(const BamRecordImpl& recordImpl) { Write(BamRecord(recordImpl)); }

    bool calculateBins_;
    std::unique_ptr<samFile, HtslibFileDeleter> file_;
    std::shared_ptr<bam_hdr_t> header_;
    std::unique_ptr<FileProducer> fileProducer_;
};

BamWriter::BamWriter(const std::string& filename, const BamHeader& header,
                     const BamWriter::CompressionLevel compressionLevel, const size_t numThreads,
                     const BinCalculationMode binCalculationMode, const bool useTempFile)
    : IRecordWriter()
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_ = std::make_unique<BamWriterPrivate>(filename, BamHeaderMemory::MakeRawHeader(header),
                                            compressionLevel, numThreads, binCalculationMode,
                                            useTempFile);
}

BamWriter::BamWriter(const std::string& filename, const BamHeader& header,
                     const BamWriter::Config& config)
    : BamWriter{filename,
                header,
                config.compressionLevel,
                config.numThreads,
                config.binCalculationMode,
                config.useTempFile}
{
}

BamWriter::BamWriter(BamWriter&&) = default;

BamWriter& BamWriter::operator=(BamWriter&&) = default;

BamWriter::~BamWriter()
{
    const auto ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    UNUSED(ret);
}

void BamWriter::TryFlush()
{
    // TODO: sanity checks on file_ & fp
    const auto ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    if (ret != 0) throw std::runtime_error{"BamWriter: could not flush output buffer contents"};
}

void BamWriter::Write(const BamRecord& record) { d_->Write(record); }

void BamWriter::Write(const BamRecord& record, int64_t* vOffset) { d_->Write(record, vOffset); }

void BamWriter::Write(const BamRecordImpl& recordImpl) { d_->Write(recordImpl); }

}  // namespace BAM
}  // namespace PacBio

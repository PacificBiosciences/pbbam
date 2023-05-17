#include "PbbamInternalConfig.h"

#include <pbbam/BamWriter.h>

#include <pbbam/BamFile.h>
#include <pbbam/Deleters.h>
#include <pbbam/Validator.h>
#include "Autovalidate.h"
#include "ErrnoReason.h"
#include "FileProducer.h"
#include "MemoryUtils.h"

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>

#include <cassert>
#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace BAM {

struct BamWriterException : public std::exception
{
    BamWriterException(std::string filename, std::string reason) : std::exception{}
    {
        std::ostringstream s;
        s << "[pbbam] BAM writer ERROR: " << reason << ":\n"
          << "  file: " << filename;
        MaybePrintErrnoReason(s);
        msg_ = s.str();
    }

    const char* what() const noexcept override { return msg_.c_str(); }

    std::string msg_;
};

class BamWriter::BamWriterPrivate
{
public:
    BamWriterPrivate(const std::string& filename, const std::shared_ptr<bam_hdr_t> rawHeader,
                     const BamWriter::CompressionLevel compressionLevel,
                     const std::size_t numThreads,
                     const BamWriter::BinCalculationMode binCalculationMode, const bool useTempFile)
        : calculateBins_{binCalculationMode == BamWriter::BinCalculation_ON}, header_{rawHeader}
    {
        if (!header_) {
            throw BamWriterException{filename, "null header provided"};
        }

        if (useTempFile) {
            fileProducer_ = std::make_unique<FileProducer>(filename);
        }

        // open file
        outputFilename_ = (fileProducer_ ? fileProducer_->TempFilename() : filename);
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel));
        file_.reset(sam_open(outputFilename_.c_str(), mode.c_str()));
        if (!file_) {
            throw BamWriterException{outputFilename_, "could not open file for writing"};
        }

        // if no explicit thread count given, attempt built-in check
        std::size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) {
                actualNumThreads = 1;
            }
        }

        // if multithreading requested, enable it
        if (actualNumThreads > 1) {
            hts_set_threads(file_.get(), actualNumThreads);
        }

        // write header
        const auto ret = sam_hdr_write(file_.get(), header_.get());
        if (ret != 0) {
            throw BamWriterException{outputFilename_, "could not write header"};
        }
    }

    void Write(const BamRecord& record)
    {
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif

        const auto& rawRecord = BamRecordMemory::GetRawData(record);

        // (probably) store bins
        // min_shift=14 & n_lvls=5 are BAM "magic numbers"
        if (calculateBins_) {
            rawRecord->core.bin =
                hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);
        }

        // write record to file
        const auto ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) {
            throw BamWriterException{outputFilename_, "could not write record"};
        }
    }

    void Write(const BamRecord& record, int64_t* vOffset)
    {
        BGZF* bgzf = file_.get()->fp.bgzf;
        assert(bgzf);
        assert(vOffset);

        // ensure offsets up-to-date
        const auto ret = bgzf_flush(bgzf);
        std::ignore = ret;

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
    std::string outputFilename_;
};

BamWriter::BamWriter(const std::string& filename, const BamHeader& header,
                     const BamWriter::CompressionLevel compressionLevel,
                     const std::size_t numThreads, const BinCalculationMode binCalculationMode,
                     const bool useTempFile)
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
{}

BamWriter::BamWriter(BamWriter&&) noexcept = default;

BamWriter& BamWriter::operator=(BamWriter&&) noexcept = default;

BamWriter::~BamWriter()
{
    const auto ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    std::ignore = ret;
}

void BamWriter::TryFlush()
{
    const auto ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    if (ret != 0) {
        throw BamWriterException{d_->outputFilename_, "could not flush buffer contents"};
    }
}

void BamWriter::Write(const BamRecord& record) { d_->Write(record); }

void BamWriter::Write(const BamRecord& record, int64_t* vOffset) { d_->Write(record, vOffset); }

void BamWriter::Write(const BamRecordImpl& recordImpl) { d_->Write(recordImpl); }

}  // namespace BAM
}  // namespace PacBio

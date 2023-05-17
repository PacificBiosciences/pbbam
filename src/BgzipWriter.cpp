#include "PbbamInternalConfig.h"

#include <pbbam/BgzipWriter.h>

#include <pbbam/Deleters.h>
#include "ErrnoReason.h"
#include "FileProducer.h"

#include <htslib/bgzf.h>
#include <htslib/hts.h>

#include <sstream>
#include <stdexcept>
#include <thread>
#include <type_traits>

namespace PacBio {
namespace BAM {

class BgzipWriter::BgzipWriterPrivate
{
public:
    BgzipWriterPrivate(std::string filename,
                       const BgzipWriterConfig& config = BgzipWriterConfig{})  // : filename_
    {
        if (config.UseTempFile) {
            fileProducer_ = std::make_unique<FileProducer>(filename);
        }

        // open file
        usingFilename_ = (fileProducer_ ? fileProducer_->TempFilename() : filename);
        const auto mode = "wb" + std::to_string(config.CompressionLevel);
        bgzf_.reset(bgzf_open(usingFilename_.c_str(), mode.c_str()));
        if (!bgzf_) {
            std::ostringstream s;
            s << "[pbbam] bgzipped file writer ERROR: could not open file for writing:\n"
              << "  file: " << usingFilename_;
            MaybePrintErrnoReason(s);
            throw std::runtime_error{s.str()};
        }

        // if no explicit thread count given, attempt built-in check
        std::size_t actualNumThreads = config.NumThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) {
                actualNumThreads = 1;
            }
        }

        // if multithreading requested, enable it
        if (actualNumThreads > 1) {
            bgzf_mt(bgzf_.get(), actualNumThreads, 256);
        }
    }

    std::size_t Write(const void* data, std::size_t numBytes)
    {
        const auto written = bgzf_write(bgzf_.get(), data, numBytes);
        if (written < 0) {
            std::ostringstream s;
            s << "[pbbam] bgzipped file writer ERROR: failed writing:\n"
              << "  file: " << usingFilename_;
            MaybePrintErrnoReason(s);
            throw std::runtime_error{s.str()};
        }
        return static_cast<std::size_t>(written);
    }

    std::string usingFilename_;
    std::unique_ptr<FileProducer> fileProducer_;
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf_;
};

BgzipWriter::BgzipWriter(std::string filename)
    : d_{std::make_unique<BgzipWriterPrivate>(std::move(filename))}
{}

BgzipWriter::BgzipWriter(std::string filename, const BgzipWriterConfig& config)
    : d_{std::make_unique<BgzipWriterPrivate>(std::move(filename), config)}
{}

BgzipWriter::BgzipWriter(BgzipWriter&&) noexcept = default;

BgzipWriter& BgzipWriter::operator=(BgzipWriter&&) noexcept = default;

BgzipWriter::~BgzipWriter() = default;

size_t BgzipWriter::Write(const void* data, std::size_t numBytes)
{
    return d_->Write(data, numBytes);
}

size_t BgzipWriter::Write(const std::string& data) { return d_->Write(data.c_str(), data.size()); }

}  // namespace BAM
}  // namespace PacBio

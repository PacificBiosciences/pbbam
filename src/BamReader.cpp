#include "PbbamInternalConfig.h"

#include <pbbam/BamReader.h>

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <sstream>
#include <stdexcept>
#include <string>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <boost/optional.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/Deleters.h>
#include <pbbam/Validator.h>

#include "Autovalidate.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

class BamReader::BamReaderPrivate
{
public:
    explicit BamReaderPrivate(std::string fn) : filename_{std::move(fn)}
    {
        auto displayFilename = [&]() {
            if (filename_ == "-")
                return std::string{" stdin"};
            else
                return "\n  file: " + filename_;
        };

        htsFile_.reset(sam_open(filename_.c_str(), "rb"));
        if (!htsFile_ || !htsFile_->fp.bgzf) {
            std::ostringstream s;
            s << "[pbbam] BAM reader ERROR: could not open for reading:" << displayFilename();
            throw std::runtime_error{s.str()};
        }

        // Use environment variable as number of BamReader threads
        if (const char* envCStr = std::getenv("PB_BAMREADER_THREADS")) {
            try {
                const int32_t numThreads = std::stoi(envCStr);
                if (numThreads <= 0) {
                    std::ostringstream s;
                    s << "[pbbam] BAM reader ERROR: environment variable PB_BAMREADER_THREADS is "
                         "not a positive, non-negative number:"
                      << envCStr;
                    throw std::runtime_error{s.str()};
                }
                p.pool = hts_tpool_init(numThreads);
                hts_set_opt(htsFile_.get(), HTS_OPT_THREAD_POOL, &p);
            } catch (const std::exception&) {
                std::ostringstream s;
                s << "[pbbam] BAM reader ERROR: environment variable PB_BAMREADER_THREADS is not a "
                     "number:"
                  << envCStr;
                throw std::runtime_error{s.str()};
            }
        }

        const auto bgzfPos = bgzf_tell(htsFile_->fp.bgzf);
        if (bgzfPos != 0) {
            std::ostringstream s;
            s << "[pbbam] BAM reader ERROR: could not read from empty input:" << displayFilename();
            throw std::runtime_error{s.str()};
        }

        const std::unique_ptr<bam_hdr_t, HtslibHeaderDeleter> hdr(sam_hdr_read(htsFile_.get()));
        if (!hdr) {
            std::ostringstream s;
            s << "[pbbam] BAM reader ERROR: could not read header from:" << displayFilename();
            throw std::runtime_error{s.str()};
        }

        header_ = BamHeaderMemory::FromRawData(hdr.get());
    }

    ~BamReaderPrivate()
    {
        if (p.pool) hts_tpool_destroy(p.pool);
    }

    htsThreadPool p = {NULL, 0};
    std::string filename_;
    std::unique_ptr<samFile, HtslibFileDeleter> htsFile_;
    BamHeader header_;
};

BamReader::BamReader() : internal::IQuery{}, d_{std::make_unique<BamReaderPrivate>("-")} {}

BamReader::BamReader(std::string fn)
    : internal::IQuery{}, d_{std::make_unique<BamReaderPrivate>(std::move(fn))}
{
}

BamReader::BamReader(BamFile bamFile) : BamReader{bamFile.Filename()} {}

BamReader::~BamReader() = default;

BGZF* BamReader::Bgzf() const { return d_->htsFile_->fp.bgzf; }

const std::string& BamReader::Filename() const { return d_->filename_; }

const BamHeader& BamReader::Header() const { return d_->header_; }

bool BamReader::GetNext(BamRecord& record)
{
    assert(BamRecordMemory::GetRawData(record).get());

    const auto result = ReadRawData(d_->htsFile_.get(), BamRecordMemory::GetRawData(record).get());

    // success
    if (result >= 0) {
        BamRecordMemory::UpdateRecordTags(record);
        record.header_ = d_->header_;
        record.ResetCachedPositions();

#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        return true;
    }

    // EOF or end-of-data range (not an error)
    else if (result == -1)
        return false;

    // error corrupted file
    else {
        std::ostringstream msg;
        msg << "[pbbam] BAM reader ERROR: cannot read from corrupted file:\n"
            << "  file: " << Filename() << '\n'
            << "  reason: ";
        if (result == -2)
            msg << "probably truncated";
        else if (result == -3)
            msg << "could not read BAM record's' core data";
        else if (result == -4)
            msg << "could not read BAM record's' variable-length data";
        else
            msg << "unknown reason (status code = " << result << ") (" << Filename() << ')';
        throw std::runtime_error{msg.str()};
    }
}

int BamReader::ReadRawData(samFile* file, bam1_t* b) { return bam_read1(file->fp.bgzf, b); }

void BamReader::VirtualSeek(int64_t virtualOffset)
{
    const auto result = bgzf_seek(Bgzf(), virtualOffset, SEEK_SET);
    if (result != 0) {
        std::ostringstream msg;
        msg << "[pbbam] BAM reader ERROR: failed to seek:\n"
            << "  file: " << Filename() << '\n'
            << "  vOffset: " << virtualOffset;
        throw std::runtime_error{msg.str()};
    }
}

int64_t BamReader::VirtualTell() const { return bgzf_tell(Bgzf()); }

}  // namespace BAM
}  // namespace PacBio

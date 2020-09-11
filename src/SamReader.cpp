#include "PbbamInternalConfig.h"

#include <pbbam/SamReader.h>

#include <cassert>
#include <cstdint>
#include <cstdio>

#include <sstream>
#include <stdexcept>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <pbbam/Deleters.h>

#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

class SamReader::SamReaderPrivate
{
public:
    explicit SamReaderPrivate(std::string fn) : filename_{std::move(fn)}
    {
        auto displayFilename = [&]() {
            if (filename_ == "-")
                return std::string{" stdin"};
            else
                return "\n  file: " + filename_;
        };

        htsFile_.reset(sam_open(filename_.c_str(), "r"));
        if (!htsFile_ || !htsFile_->fp.hfile) {
            std::ostringstream s;
            s << "[pbbam] SAM reader ERROR: could not open:" << displayFilename();
            throw std::runtime_error{s.str()};
        }

        char c;
        const auto peek = hpeek(htsFile_->fp.hfile, &c, 1);
        if (peek == 0) {
            std::ostringstream s;
            s << "[pbbam] SAM reader ERROR: could not read from empty input:" << displayFilename();
            throw std::runtime_error{s.str()};
        }

        hdr_.reset(sam_hdr_read(htsFile_.get()));
        if (!hdr_ || hdr_->l_text == 0) {
            std::ostringstream s;
            s << "[pbbam] SAM reader ERROR: could not read header from:" << displayFilename();
            throw std::runtime_error{s.str()};
        }
        fullHeader_ = BamHeaderMemory::FromRawData(hdr_.get());
    }

    std::string filename_;
    std::unique_ptr<samFile, HtslibFileDeleter> htsFile_;
    std::unique_ptr<bam_hdr_t, HtslibHeaderDeleter> hdr_;
    BamHeader fullHeader_;  // fully parsed, for `GetNext`
};

SamReader::SamReader() : internal::IQuery{}, d_{std::make_unique<SamReaderPrivate>("-")} {}

SamReader::SamReader(std::string fn)
    : internal::IQuery{}, d_{std::make_unique<SamReaderPrivate>(std::move(fn))}
{
}

SamReader::~SamReader() = default;

const std::string& SamReader::Filename() const { return d_->filename_; }

const BamHeader& SamReader::Header() const { return d_->fullHeader_; }

bool SamReader::GetNext(BamRecord& record)
{
    auto* b = BamRecordMemory::GetRawData(record).get();
    assert(b);
    const auto result = sam_read1(d_->htsFile_.get(), d_->hdr_.get(), b);

    // success
    if (result >= 0) {
        BamRecordMemory::UpdateRecordTags(record);
        record.header_ = d_->fullHeader_;
        record.ResetCachedPositions();
        return true;
    }

    // EOF or end-of-data range (not an error)
    else if (result == -1)
        return false;

    // error corrupted file
    else {
        std::ostringstream msg;
        msg << "[pbbam] SAM reader ERROR: cannot read from corrupted file:\n"
            << "  file: " << Filename() << '\n'
            << "  reason: ";
        if (result == -2)
            msg << "probably truncated";
        else if (result == -3)
            msg << "could not read SAM record's' core data";
        else if (result == -4)
            msg << "could not read SAM record's' variable-length data";
        else
            msg << "unknown reason (status code = " << result << ") (" << Filename() << ')';
        throw std::runtime_error{msg.str()};
    }
}

}  // namespace BAM
}  // namespace PacBio

#include "PbbamInternalConfig.h"

#include <pbbam/BamFile.h>

#include <sys/stat.h>

#include <cassert>
#include <cstdint>

#include <memory>
#include <sstream>
#include <string>

#include <htslib/sam.h>

#include <pbbam/Deleters.h>
#include <pbbam/PbiFile.h>

#include "Autovalidate.h"
#include "ErrnoReason.h"
#include "FileUtils.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

class BamFile::BamFilePrivate
{
public:
    explicit BamFilePrivate(std::string fn) : filename_{std::move(fn)}
    {
        // attempt open
        auto f = RawOpen();

#if !defined(PBBAM_NO_CHECK_EOF) || PBBAM_AUTOVALIDATE
        // sanity check on file
        const auto eofCheck = bgzf_check_EOF(f->fp.bgzf);
        if (eofCheck <= 0) {
            // 1:  EOF present & correct
            // 2:  not seekable (e.g. reading from stdin)
            // 0:  EOF absent
            // -1: some other error
            std::ostringstream e;
            if (eofCheck == 0) {
                e << "[pbbam] BAM file ERROR: missing EOF block:\n"
                  << "  file: " << fn;
            } else {
                e << "[pbbam] BAM file ERROR: unknown error encountered while checking EOF:\n"
                  << "  file: " << fn;
                MaybePrintErrnoReason(e);
                e << "\n  htslib status code: " << eofCheck;
                throw std::runtime_error{e.str()};
            }
        }
#endif

        // attempt fetch header
        std::unique_ptr<bam_hdr_t, HtslibHeaderDeleter> hdr(sam_hdr_read(f.get()));
        header_ = BamHeaderMemory::FromRawData(hdr.get());
    }

    std::unique_ptr<BamFilePrivate> DeepCopy()
    {
        return std::make_unique<BamFilePrivate>(filename_);
    }

    bool HasEOF() const
    {
        // streamed input is unknown, since it's not random-accessible
        if (filename_ == "-") {
            return false;
        }

        // attempt open
        const auto f = RawOpen();
        return RawEOFCheck(f) == 1;
    }

    int RawEOFCheck(const std::unique_ptr<samFile, HtslibFileDeleter>& f) const
    {
        assert(f);
        assert(f->fp.bgzf);
        return bgzf_check_EOF(f->fp.bgzf);
    }

    std::unique_ptr<samFile, HtslibFileDeleter> RawOpen() const
    {
        std::unique_ptr<samFile, HtslibFileDeleter> f(sam_open(filename_.c_str(), "rb"));
        if (!f || !f->fp.bgzf) {
            std::ostringstream s;
            s << "[pbbam] BAM file ERROR: could not open:\n"
              << "  file: " << filename_;
            MaybePrintErrnoReason(s);
            throw std::runtime_error{s.str()};
        }
        if (f->format.format != bam) {
            std::ostringstream s;
            s << "[pbbam] BAM file ERROR: expected BAM, encountered different format:\n"
              << "  file: " << filename_;
            throw std::runtime_error{s.str()};
        }

        return f;
    }

    std::string filename_;
    BamHeader header_;
    int64_t firstAlignmentOffset_;
};

BamFile::BamFile(std::string filename) : d_{std::make_unique<BamFilePrivate>(std::move(filename))}
{
}

BamFile::BamFile(const BamFile& other) : d_{other.d_->DeepCopy()} {}

BamFile::BamFile(BamFile&&) noexcept = default;

BamFile& BamFile::operator=(const BamFile& other)
{
    if (this != &other) {
        d_ = other.d_->DeepCopy();
    }
    return *this;
}

BamFile& BamFile::operator=(BamFile&&) noexcept = default;

BamFile::~BamFile() = default;

void BamFile::CreatePacBioIndex() const { PbiFile::CreateFrom(*this); }

void BamFile::CreateStandardIndex() const
{
    const auto ret = bam_index_build(d_->filename_.c_str(), 0);
    if (ret != 0) {
        std::ostringstream s;
        s << "[pbbam] BAM file ERROR: could not create *.bai index:\n"
          << "  file: " << d_->filename_;
        MaybePrintErrnoReason(s);
        s << "\n  htslib status code: " << ret;
        throw std::runtime_error{s.str()};
    }
}

void BamFile::EnsurePacBioIndexExists() const
{
    if (!PacBioIndexExists()) {
        CreatePacBioIndex();
    }
}

void BamFile::EnsureStandardIndexExists() const
{
    if (!StandardIndexExists()) {
        CreateStandardIndex();
    }
}

const std::string& BamFile::Filename() const { return d_->filename_; }

bool BamFile::HasEOF() const { return d_->HasEOF(); }

bool BamFile::HasReference(const std::string& name) const { return d_->header_.HasSequence(name); }

const BamHeader& BamFile::Header() const { return d_->header_; }

bool BamFile::IsPacBioBAM() const { return !d_->header_.PacBioBamVersion().empty(); }

bool BamFile::PacBioIndexExists() const { return FileUtils::Exists(PacBioIndexFilename()); }

std::string BamFile::PacBioIndexFilename() const { return d_->filename_ + ".pbi"; }

bool BamFile::PacBioIndexIsNewer() const
{
    const auto bamTimestamp = FileUtils::LastModified(Filename());
    const auto pbiTimestamp = FileUtils::LastModified(PacBioIndexFilename());
    return bamTimestamp <= pbiTimestamp;
}

int BamFile::ReferenceId(const std::string& name) const { return d_->header_.SequenceId(name); }

uint32_t BamFile::ReferenceLength(const std::string& name) const
{
    return ReferenceLength(ReferenceId(name));
}

uint32_t BamFile::ReferenceLength(const int id) const
{
    return std::stoul(d_->header_.SequenceLength(id));
}

std::string BamFile::ReferenceName(const int id) const { return d_->header_.SequenceName(id); }

bool BamFile::StandardIndexExists() const { return FileUtils::Exists(StandardIndexFilename()); }

std::string BamFile::StandardIndexFilename() const { return d_->filename_ + ".bai"; }

bool BamFile::StandardIndexIsNewer() const
{
    const auto bamTimestamp = FileUtils::LastModified(Filename());
    const auto baiTimestamp = FileUtils::LastModified(StandardIndexFilename());
    return bamTimestamp <= baiTimestamp;
}

}  // namespace BAM
}  // namespace PacBio

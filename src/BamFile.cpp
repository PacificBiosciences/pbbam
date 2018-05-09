// File Description
/// \file BamFile.cpp
/// \brief Implements the BamFile class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/BamFile.h"

#include <sys/stat.h>
#include <cassert>
#include <cstdint>
#include <memory>
#include <sstream>

#include <htslib/sam.h>

#include "Autovalidate.h"
#include "FileUtils.h"
#include "MemoryUtils.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiFile.h"

namespace PacBio {
namespace BAM {
namespace internal {

class BamFilePrivate
{
public:
    BamFilePrivate(std::string fn) : filename_{std::move(fn)}, firstAlignmentOffset_{-1}
    {
        // ensure we've updated htslib verbosity with requested verbosity here
        hts_verbose = (PacBio::BAM::HtslibVerbosity == -1 ? 0 : PacBio::BAM::HtslibVerbosity);

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
            if (eofCheck == 0)
                e << fn << " : is missing EOF block\n";
            else
                e << fn << " : unknown error while checking EOF block\n";
            throw std::runtime_error{e.str()};
        }
#endif

        // attempt fetch header
        std::unique_ptr<bam_hdr_t, internal::HtslibHeaderDeleter> hdr(sam_hdr_read(f.get()));
        header_ = internal::BamHeaderMemory::FromRawData(hdr.get());

        // cache first alignment offset
        firstAlignmentOffset_ = bgzf_tell(f->fp.bgzf);
    }

    std::unique_ptr<BamFilePrivate> DeepCopy()
    {
        return std::make_unique<BamFilePrivate>(filename_);
    }

    bool HasEOF() const
    {
        // streamed input is unknown, since it's not random-accessible
        if (filename_ == "-") return false;

        // attempt open
        auto f = RawOpen();
        return RawEOFCheck(f) == 1;
    }

    int RawEOFCheck(const std::unique_ptr<samFile, internal::HtslibFileDeleter>& f) const
    {
        assert(f);
        assert(f->fp.bgzf);
        return bgzf_check_EOF(f->fp.bgzf);
    }

    std::unique_ptr<samFile, internal::HtslibFileDeleter> RawOpen() const
    {
        std::unique_ptr<samFile, internal::HtslibFileDeleter> f(sam_open(filename_.c_str(), "rb"));
        if (!f || !f->fp.bgzf) throw std::runtime_error{"could not open BAM file: " + filename_};
        if (f->format.format != bam) throw std::runtime_error{"expected BAM, unknown format"};
        return f;
    }

public:
    std::string filename_;
    BamHeader header_;
    int64_t firstAlignmentOffset_;
};

}  // namespace internal

// ------------------------
// BamFile implementation
// ------------------------

BamFile::BamFile(std::string filename)
    : d_{std::make_unique<internal::BamFilePrivate>(std::move(filename))}
{
}

BamFile::BamFile(const BamFile& other) : d_{other.d_->DeepCopy()} {}

BamFile::BamFile(BamFile&& other) : d_{std::move(other.d_)} {}

BamFile& BamFile::operator=(const BamFile& other)
{
    if (this != &other) {
        d_ = other.d_->DeepCopy();
    }
    return *this;
}

BamFile& BamFile::operator=(BamFile&& other)
{
    if (this != &other) {
        d_ = std::move(other.d_);
    }
    return *this;
}

BamFile::~BamFile() {}

void BamFile::CreatePacBioIndex() const { PbiFile::CreateFrom(*this); }

void BamFile::CreateStandardIndex() const
{
    if (bam_index_build(d_->filename_.c_str(), 0) != 0)
        throw std::runtime_error{"could not build BAI index"};
}

void BamFile::EnsurePacBioIndexExists() const
{
    if (!PacBioIndexExists()) CreatePacBioIndex();
}

void BamFile::EnsureStandardIndexExists() const
{
    if (!StandardIndexExists()) CreateStandardIndex();
}

const std::string& BamFile::Filename() const { return d_->filename_; }

int64_t BamFile::FirstAlignmentOffset() const { return d_->firstAlignmentOffset_; }

bool BamFile::HasEOF() const { return d_->HasEOF(); }

bool BamFile::HasReference(const std::string& name) const { return d_->header_.HasSequence(name); }

const BamHeader& BamFile::Header() const { return d_->header_; }

bool BamFile::IsPacBioBAM() const { return !d_->header_.PacBioBamVersion().empty(); }

bool BamFile::PacBioIndexExists() const
{
    return internal::FileUtils::Exists(PacBioIndexFilename());
}

std::string BamFile::PacBioIndexFilename() const { return d_->filename_ + ".pbi"; }

bool BamFile::PacBioIndexIsNewer() const
{
    const auto bamTimestamp = internal::FileUtils::LastModified(Filename());
    const auto pbiTimestamp = internal::FileUtils::LastModified(PacBioIndexFilename());
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

bool BamFile::StandardIndexExists() const
{
    return internal::FileUtils::Exists(StandardIndexFilename());
}

std::string BamFile::StandardIndexFilename() const { return d_->filename_ + ".bai"; }

bool BamFile::StandardIndexIsNewer() const
{
    const auto bamTimestamp = internal::FileUtils::LastModified(Filename());
    const auto baiTimestamp = internal::FileUtils::LastModified(StandardIndexFilename());
    return bamTimestamp <= baiTimestamp;
}

}  // namespace BAM
}  // namespace PacBio

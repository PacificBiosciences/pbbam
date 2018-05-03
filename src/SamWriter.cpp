// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SamWriter.h"

#include <htslib/hfile.h>
#include <htslib/sam.h>

#include "Autovalidate.h"
#include "FileProducer.h"
#include "MemoryUtils.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/Validator.h"

namespace PacBio {
namespace BAM {
namespace internal {

class SamWriterPrivate : public internal::FileProducer
{
public:
    SamWriterPrivate(std::string filename, const std::shared_ptr<bam_hdr_t> rawHeader)
        : internal::FileProducer{std::move(filename)}, header_{rawHeader}
    {
        if (!header_) throw std::runtime_error{"null header"};

        // open file
        const auto& usingFilename = TempFilename();
        const std::string mode(1, 'w');
        file_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!file_)
            throw std::runtime_error{"could not open file: " + usingFilename + "for writing"};

        // write header
        const auto ret = sam_hdr_write(file_.get(), header_.get());
        if (ret != 0) throw std::runtime_error{"could not write header"};
    }

    void TryFlush();
    void Write(const BamRecord& record);

private:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> file_;
    std::shared_ptr<bam_hdr_t> header_;
};

void SamWriterPrivate::TryFlush()
{
    const auto ret = file_.get()->fp.hfile;
    if (ret != nullptr) throw std::runtime_error{"could not flush output buffer contents"};
}

void SamWriterPrivate::Write(const BamRecord& record)
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(record);
#endif

    const auto rawRecord = internal::BamRecordMemory::GetRawData(record);

    // store bin number
    // min_shift=14 & n_lvls=5 are SAM/BAM "magic numbers"
    rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

    // write record to file
    const int ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
    if (ret <= 0) throw std::runtime_error("could not write record");
}

}  // namespace internal

SamWriter::SamWriter(std::string filename, const BamHeader& header)
    : IRecordWriter()
    , d_{std::make_unique<internal::SamWriterPrivate>(
          std::move(filename), internal::BamHeaderMemory::MakeRawHeader(header))}
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
}

SamWriter::~SamWriter() {}

void SamWriter::TryFlush() { d_->TryFlush(); }

void SamWriter::Write(const BamRecord& record) { d_->Write(record); }

void SamWriter::Write(const BamRecordImpl& recordImpl) { d_->Write(BamRecord{recordImpl}); }

}  // namespace BAM
}  // namespace PacBio

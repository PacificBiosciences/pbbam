// File Description
/// \file IndexedBamWriter.cpp
/// \brief Implements the IndexedBamWriter class
//
// Author: Derek Barnett

#include "pbbam/IndexedBamWriter.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include "pbbam/BamHeader.h"
#include "pbbam/BamRecord.h"
#include "pbbam/MakeUnique.h"
#include "pbbam/PbiBuilder.h"
#include "pbbam/Unused.h"
#include "pbbam/Validator.h"

#include "FileProducer.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {
namespace internal {

class IndexedBamWriterPrivate : public internal::FileProducer
{
public:
    IndexedBamWriterPrivate(const std::string& outputFilename, std::shared_ptr<bam_hdr_t> rawHeader)
        : internal::FileProducer{outputFilename}
        , header_{rawHeader}
        , builder_{outputFilename + ".pbi"}
        , previousBlockAddress_{0}
    {
        if (!header_) throw std::runtime_error{"null header"};

        // open file
        const auto& usingFilename = TempFilename();
        file_.reset(sam_open(usingFilename.c_str(), "wb"));
        if (!file_)
            throw std::runtime_error{"could not open file" + usingFilename + " for writing"};

        // write header
        const auto ret = sam_hdr_write(file_.get(), header_.get());
        if (ret != 0) throw std::runtime_error{"could not write header"};

        // store first alignment block
        previousBlockAddress_ = file_.get()->fp.bgzf->block_address;
    }

    ~IndexedBamWriterPrivate()
    {
        // ensure last remaining bits are flushed to file
        const auto ret = bgzf_flush(file_.get()->fp.bgzf);
        UNUSED(ret);
    }

public:
    void Write(const BamRecord& record)
    {
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        const auto rawRecord = internal::BamRecordMemory::GetRawData(record);
        BGZF* bgzf = file_.get()->fp.bgzf;
        assert(bgzf);

        //
        // Fetch record's start offset.
        //
        // If we're still in the same block from the last record written, we
        // need to flush to get the proper offset.
        //
        if (bgzf->block_address == previousBlockAddress_) {
            const auto ret = bgzf_flush(bgzf);
            UNUSED(ret);
        }
        const int64_t vOffset = bgzf_tell(bgzf);

        // update bin
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file & PBI builder
        const auto ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw std::runtime_error{"could not write record"};
        builder_.AddRecord(record, vOffset);

        // update block address
        previousBlockAddress_ = bgzf->block_address;
    }

public:
    std::unique_ptr<samFile, internal::HtslibFileDeleter> file_;
    std::shared_ptr<bam_hdr_t> header_;
    PbiBuilder builder_;
    int64_t previousBlockAddress_;
};

}  // namespace internal

IndexedBamWriter::IndexedBamWriter(const std::string& outputFilename, const BamHeader& header)
    : IRecordWriter()
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_ = std::make_unique<internal::IndexedBamWriterPrivate>(
        outputFilename, internal::BamHeaderMemory::MakeRawHeader(header));
}

IndexedBamWriter::~IndexedBamWriter() {}

void IndexedBamWriter::TryFlush()
{
    const auto ret = bgzf_flush(d_->file_.get()->fp.bgzf);
    if (ret != 0) throw std::runtime_error{"could not flush output buffer contents"};
}

void IndexedBamWriter::Write(const BamRecord& record) { d_->Write(record); }

void IndexedBamWriter::Write(const BamRecordImpl& record) { d_->Write(BamRecord{record}); }

}  // namespace BAM
}  // namespace PacBio

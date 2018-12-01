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

class SamWriter::SamWriterPrivate : public FileProducer
{
public:
    SamWriterPrivate(std::string filename, const std::shared_ptr<bam_hdr_t> rawHeader)
        : FileProducer{std::move(filename)}, header_{rawHeader}
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

    void Write(BamRecord record)
    {
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif

        const auto rawRecord = BamRecordMemory::GetRawData(record);

        // store bin number
        // min_shift=14 & n_lvls=5 are SAM/BAM "magic numbers"
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // Maybe adjust location of long CIGAR (>65535 ops) data, depending on the
        // runtime htslib version.
        //
        // SAM formatting in htslib verions previous to 1.7 are unaware of the new
        // long CIGAR implementation ("CG") tag. So we need to move that back to the
        // "standard" field so that SAM output is correct. Versions >=1.7 properly
        // display long CIGARs.
        //
        // This transform will become unecessary when we drop support for htslib pre-v1.7.
        //
        static const bool has_native_long_cigar_support = DoesHtslibSupportLongCigar();
        const auto cigar = record.CigarData();
        if (!has_native_long_cigar_support && cigar.size() > 65535) {
            if (record.Impl().HasTag("CG")) record.Impl().RemoveTag("CG");
            record.Impl().SetCigarData(cigar);
        }

        // write record to file
        const int ret = sam_write1(file_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw std::runtime_error("could not write record");
    }

    std::unique_ptr<samFile, HtslibFileDeleter> file_;
    std::shared_ptr<bam_hdr_t> header_;
};

SamWriter::SamWriter(std::string filename, const BamHeader& header)
    : IRecordWriter()
    , d_{std::make_unique<SamWriterPrivate>(std::move(filename),
                                            BamHeaderMemory::MakeRawHeader(header))}
{
#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
}

SamWriter::~SamWriter() {}

void SamWriter::TryFlush()
{
    const auto ret = d_->file_.get()->fp.hfile;
    if (ret != nullptr) throw std::runtime_error{"could not flush output buffer contents"};
}

void SamWriter::Write(const BamRecord& record) { d_->Write(record); }

void SamWriter::Write(const BamRecordImpl& recordImpl) { Write(BamRecord{recordImpl}); }

}  // namespace BAM
}  // namespace PacBio

#include "PbbamInternalConfig.h"

#include <pbbam/IndexedBamWriter.h>

#include <pbbam/BamHeader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordImpl.h>
#include <pbbam/BamWriter.h>
#include <pbbam/Deleters.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/RecordType.h>
#include <pbbam/Validator.h>
#include "ErrnoReason.h"
#include "FileProducer.h"
#include "MemoryUtils.h"
#include "PbiBuilderBase.h"

#include <pbcopper/utility/Deleters.h>

#include <boost/numeric/conversion/cast.hpp>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

#include <array>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <type_traits>

#include <cstdint>

#include <sys/stat.h>

namespace PacBio {
namespace BAM {

// using PbiBuilderException = PbiBuilderException;
// using IndexedBamWriterException = IndexedBamWriterException;

struct GzIndexEntry
{
    int64_t vAddress;
    int64_t uAddress;
};

// TODO: come back to refseqs, sorting, etc
class PbiBuilder2 : public PacBio::BAM::PbiBuilderBase
{
public:
    PbiBuilder2(const std::string& bamFilename, const std::string& pbiFilename,
                const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads,
                const size_t fileBufferSize)
        //                const size_t numReferenceSequences = 0
        //                const bool isCoordinateSorted = false
        : PacBio::BAM::PbiBuilderBase{pbiFilename, compressionLevel, numThreads, fileBufferSize}
        , bamFilename_{bamFilename}
    {}

    std::vector<GzIndexEntry> LoadGzi()
    {
        //
        // Open GZI file & load its contents. About to use for offset transformation.
        //

        const std::string gziFn{bamFilename_ + ".gzi"};
        std::unique_ptr<FILE, Utility::FileDeleter> gziFile{fopen(gziFn.c_str(), "rb")};
        if (!gziFile) {
            throw IndexedBamWriterException{gziFn, "could not open *.gzi file"};
        }

        uint64_t numElements;
        if (fread(&numElements, sizeof(numElements), 1, gziFile.get()) < 1) {
            throw IndexedBamWriterException{gziFn, "could not read from *.gzi file"};
        }
        if (ed_is_big()) {
            ed_swap_8(numElements);
        }

        std::vector<GzIndexEntry> result;
        result.reserve(numElements);
        for (uint32_t i = 0; i < numElements; ++i) {
            GzIndexEntry entry;
            const auto vReturn = fread(&entry.vAddress, sizeof(entry.vAddress), 1, gziFile.get());
            const auto uReturn = fread(&entry.uAddress, sizeof(entry.uAddress), 1, gziFile.get());
            if (vReturn < 1 || uReturn < 1) {
                throw IndexedBamWriterException{gziFn, "could not read from *.gzi file"};
            }

            if (ed_is_big()) {
                ed_swap_8(entry.vAddress);
                ed_swap_8(entry.uAddress);
            }
            result.push_back(std::move(entry));
        }

        if (result.empty()) {
            std::ostringstream s;
            s << "[pbbam] indexed BAM writer ERROR: empty GZI index\n"
              << "  file: " << bamFilename_ + ".gzi";
            throw std::runtime_error{s.str()};
        }

        return result;
    }

    void WriteVirtualOffsets() final
    {
        auto index = LoadGzi();
        std::sort(index.begin(), index.end(),
                  [](const GzIndexEntry& lhs, const GzIndexEntry& rhs) -> bool {
                      return lhs.uAddress < rhs.uAddress;
                  });

        size_t k = 0;
        for (const auto& block : fileOffsetField_.blocks_) {
            LoadFieldBlockFromTempFile(fileOffsetField_, block);

            // transform offsets from GZI
            for (size_t j = 0; j < fileOffsetField_.buffer_.size(); ++j) {
                while ((k < index.size() - 1) && (static_cast<uint64_t>(index.at(k + 1).uAddress) <=
                                                  fileOffsetField_.buffer_[j])) {
                    ++k;
                }
                const GzIndexEntry& e = index.at(k);
                const int64_t uOffset = fileOffsetField_.buffer_[j] - e.uAddress;
                const auto result = ((e.vAddress << 16) | uOffset);
                fileOffsetField_.buffer_[j] = result;
            }
            WriteBgzfVector(pbiFile_.get(), fileOffsetField_.buffer_);
        }
    }

private:
    std::string bamFilename_;
};

// htslib >= v.10
#if defined(HTS_VERSION) && HTS_VERSION >= 101000

class IndexedBamWriter::IndexedBamWriterPrivate2  //: public internal::FileProducer
{
public:
    IndexedBamWriterPrivate2(const std::string& outputFilename, std::shared_ptr<bam_hdr_t> header,
                             const BamWriter::CompressionLevel bamCompressionLevel,
                             const size_t numBamThreads,
                             const PbiBuilder::CompressionLevel pbiCompressionLevel,
                             const size_t numPbiThreads, const size_t /*numGziThreads*/,
                             const size_t tempFileBufferSize)
        : bamFilename_{outputFilename}, header_{header}
    {
        OpenBam(bamCompressionLevel, numBamThreads);
        OpenPbi(pbiCompressionLevel, numPbiThreads, tempFileBufferSize);
        isOpen_ = true;
    }

    ~IndexedBamWriterPrivate2() noexcept
    {
        if (isOpen_) {
            try {
                Close();
            } catch (...) {
                // swallow any exceptions & remain no-throw from dtor
            }
        }
    }

    void Close()
    {
        // NOTE: keep this order of closing ( BAM -> PBI )
        CloseBam();
        ClosePbi();

        remove(std::string{bamFilename_ + ".gzi"}.c_str());
        isOpen_ = false;
    }

    void CloseBam()
    {
        auto ret = bgzf_flush(bam_.get()->fp.bgzf);

        // Dump GZI contents to disk.
        const std::string gziFn{bamFilename_ + ".gzi"};
        ret = bgzf_index_dump(bam_.get()->fp.bgzf, gziFn.c_str(), nullptr);
        std::ignore = ret;
        bam_.reset();
    }

    void ClosePbi() { builder_->Close(); }

    void OpenBam(const BamWriter::CompressionLevel compressionLevel, const size_t numThreads)
    {
        //
        // TODO: Compression level & numThreads are hardcoded here. Ok for
        //       prototyping but need to be tune-able via API.
        //

        if (!header_) {
            throw IndexedBamWriterException{bamFilename_, "null header provided"};
        }

        // open output BAM
        const auto usingFilename = bamFilename_;
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel));
        bam_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!bam_) {
            throw IndexedBamWriterException{usingFilename, "could not open file for writing"};
        }

        const auto indexInit = bgzf_index_build_init(bam_.get()->fp.bgzf);
        if (indexInit != 0) {
            throw IndexedBamWriterException{usingFilename,
                                            "could not open initialize on-the-fly gzi index"};
        }

        // maybe set multithreaded writing
        size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) {
                actualNumThreads = 1;
            }
        }
        if (actualNumThreads > 1) {
            hts_set_threads(bam_.get(), actualNumThreads);
        }

        // write header
        auto ret = sam_hdr_write(bam_.get(), header_.get());
        if (ret != 0) {
            throw IndexedBamWriterException{usingFilename, "could not write header"};
        }
        ret = bgzf_flush(bam_.get()->fp.bgzf);

        // store file positions after header
        auto headerLength = [](const bam_hdr_t* hdr) -> size_t {
            const size_t textHeader = 12 + hdr->l_text;
            size_t refHeader = 0;
            for (int i = 0; i < hdr->n_targets; ++i) {
                char* n = hdr->target_name[i];
                refHeader += (8 + (strlen(n) + 1));
            }
            return textHeader + refHeader;
        };
        uncompressedFilePos_ = headerLength(header_.get());
    }

    void OpenPbi(const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads,
                 const size_t fileBufferSize)
    {
        builder_ = std::make_unique<PbiBuilder2>(bamFilename_, bamFilename_ + ".pbi",
                                                 compressionLevel, numThreads, fileBufferSize);
    }

    void Write(const BamRecord& record)
    {
// TODO: add API to auto-skip this without special compile flag
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        // add record & its to index builder.
        //
        // NOTE: This is the record's postiion as if it were _uncompressed_. We
        //       will return with GZI data later to transform it into BAM
        //       "virtual offset".
        //
        builder_->AddRecord(record, uncompressedFilePos_);

        const auto& rawRecord = BamRecordMemory::GetRawData(record);

        // update bin
        // min_shift=14 & n_lvls=5 are BAM "magic numbers"
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file
        const auto ret = sam_write1(bam_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) {
            throw IndexedBamWriterException{bamFilename_, "could not write record"};
        }

        // update file position
        auto recordLength = [](bam1_t* b) {
            auto* c = &b->core;

            static constexpr size_t fixedLength = 36;
            const size_t qnameLength = (c->l_qname - c->l_extranul);

            // TODO: long CIGAR handling... sigh...

            size_t remainingLength = 0;
            if (c->n_cigar <= 0xffff) {
                remainingLength = (b->l_data - c->l_qname);
            } else {
                const size_t cigarEnd = ((uint8_t*)bam_get_cigar(b) - b->data) + (c->n_cigar * 4);
                remainingLength = 8 + (b->l_data - cigarEnd) + 4 + (4 * c->n_cigar);
            }

            return fixedLength + qnameLength + remainingLength;
        };
        uncompressedFilePos_ += recordLength(rawRecord.get());
    }

private:
    std::string bamFilename_;

    std::shared_ptr<bam_hdr_t> header_;
    std::unique_ptr<samFile, HtslibFileDeleter> bam_;
    std::unique_ptr<PbiBuilder2> builder_;
    bool isOpen_ = false;
    int64_t uncompressedFilePos_ = 0;
};

#else  // htslib < v1.10

class IndexedBamWriter::IndexedBamWriterPrivate2  //: public internal::FileProducer
{
public:
    IndexedBamWriterPrivate2(const std::string& outputFilename, std::shared_ptr<bam_hdr_t> header,
                             const BamWriter::CompressionLevel bamCompressionLevel,
                             const size_t numBamThreads,
                             const PbiBuilder::CompressionLevel pbiCompressionLevel,
                             const size_t numPbiThreads, const size_t numGziThreads,
                             const size_t tempFileBufferSize)
        : bamFilename_{outputFilename}, header_{header}
    {
        OpenBam(bamCompressionLevel, numBamThreads);
        OpenGzi(numGziThreads);
        OpenPbi(pbiCompressionLevel, numPbiThreads, tempFileBufferSize);
        isOpen_ = true;
    }

    ~IndexedBamWriterPrivate2() noexcept
    {
        if (isOpen_) {
            try {
                Close();
            } catch (...) {
                // swallow any exceptions & remain no-throw from dtor
            }
        }
    }

    void Close()
    {
        // NOTE: keep this order of closing ( BAM -> GZI -> PBI )
        CloseBam();
        CloseGzi();
        ClosePbi();

        remove(std::string{bamFilename_ + ".gzi"}.c_str());
        isOpen_ = false;
    }

    void CloseBam()
    {
        const auto ret = bgzf_flush(bam_.get()->fp.bgzf);
        std::ignore = ret;
        bam_.reset();
    }

    void CloseGzi()
    {
        done_ = true;
        gziThread_.join();
    }

    void ClosePbi() { builder_->Close(); }

    void OpenBam(const BamWriter::CompressionLevel compressionLevel, const size_t numThreads)
    {
        //
        // TODO: Compression level & numThreads are hardcoded here. Ok for
        //       prototyping but need to be tune-able via API.
        //

        // open output BAM
        const auto usingFilename = bamFilename_;
        const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel));
        bam_.reset(sam_open(usingFilename.c_str(), mode.c_str()));
        if (!bam_)
            throw IndexedBamWriterException{usingFilename, "could not open file for writing"};

        // maybe set multithreaded writing
        size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) actualNumThreads = 1;
        }
        if (actualNumThreads > 1) hts_set_threads(bam_.get(), actualNumThreads);

        // write header
        if (!header_) {
            std::ostringstream s;
            s << "[pbbam] indexed BAM writer ERROR: invalid header provided\n"
              << "  file: " << usingFilename;
            MaybePrintErrnoReason(s);
            throw std::runtime_error{s.str()};
        }
        auto ret = sam_hdr_write(bam_.get(), header_.get());
        if (ret != 0) throw IndexedBamWriterException{usingFilename, "could not write header"};
        ret = bgzf_flush(bam_.get()->fp.bgzf);

        // store file positions after header
        auto headerLength = [](const bam_hdr_t* hdr) -> size_t {
            const size_t textHeader = 12 + hdr->l_text;
            size_t refHeader = 0;
            for (int i = 0; i < hdr->n_targets; ++i) {
                char* n = hdr->target_name[i];
                refHeader += (8 + (strlen(n) + 1));
            }
            return textHeader + refHeader;
        };
        uncompressedFilePos_ = headerLength(header_.get());
    }

    void OpenGzi(size_t numThreads)
    {
        size_t actualNumThreads = numThreads;
        if (actualNumThreads == 0) {
            actualNumThreads = std::thread::hardware_concurrency();

            // if still unknown, default to single-threaded
            if (actualNumThreads == 0) actualNumThreads = 1;
        }
        gziThread_ = std::thread{&IndexedBamWriterPrivate2::RunGziThread, this, actualNumThreads};
    }

    void OpenPbi(const PbiBuilder::CompressionLevel compressionLevel, const size_t numThreads,
                 const size_t fileBufferSize)
    {
        builder_ = std::make_unique<PbiBuilder2>(bamFilename_, bamFilename_ + ".pbi",
                                                 compressionLevel, numThreads, fileBufferSize);
    }

    void RunGziThread(size_t numThreads)
    {
        //
        // This thread is the GZI index-enabled reader that trails the writer
        // thread(s). It checks for changes in the output BAM's file size &
        // reads whatever data is available. When writing is complete, it reads
        // anything that might remain & dumps the GZI index contents to disk.
        // This index is used downstream to generate records' "virtual offsets".
        //

        const auto& bamFilename = bamFilename_;
        std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf;

        struct stat st;
        int ret = 0;
        int64_t lastFileSize = 0;
        int64_t numBytesRead = 0;

        auto initBgzf = [&bgzf, &bamFilename, numThreads]() {
            bgzf.reset(bgzf_open(bamFilename.c_str(), "rb"));
            if (!bgzf)
                throw IndexedBamWriterException{bamFilename, "could not open trailing BAM reader"};
            bgzf_index_build_init(bgzf.get());
            if (numThreads > 1) bgzf_mt(bgzf.get(), numThreads, 256);
        };

        // main thread loop
        while (true) {
            // Quit if writer thread(s) are finished.
            if (done_) break;

            if (stat(bamFilename.c_str(), &st) != 0) {
                gziStatus_ = GziStatus::MISC_ERROR;
                return;
            }
            if (st.st_size > lastFileSize) {
                lastFileSize = st.st_size;
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                continue;
            }

            // Don't read unless we can guarantee we won't catch up to the end of the file.
            // Otherwise htslib will think the file has been truncated and throw errors.
            // This is a touch tricky because we're reading in multi-thread mode so htslib
            // will speculatively start grabbing blocks.  So we're going to stay *well* behind
            // This needs be made more robust.
            //
            // Note: It's worth noting that bgzf->block_clength might only be the length of the
            // compressed *payload*, meaning if there is any other header/metadata/etc on disk
            // in the actual file, our estimation of our trailing distance might be off.  If
            // this ever starts throwing exceptions we'll have to look more in to this...
            while (lastFileSize - numBytesRead > 100 * BGZF_MAX_BLOCK_SIZE) {
                // Open BAM reader if not already open.  Need to make sure we don't open it
                // until we've already established our trailing distance.
                if (!bgzf) initBgzf();

                auto result = bgzf_read_block(bgzf.get());
                if (result != 0) {
                    gziStatus_ = GziStatus::IO_ERROR;
                    return;
                }
                if (bgzf->block_length == 0) {
                    gziStatus_ = GziStatus::TRAIL_ERROR;
                    return;
                }
                numBytesRead += bgzf->block_clength;
            }

            // Only update if thigs have appreciably fallen behind
            if (lastFileSize - numBytesRead > 1.10 * maxTrailingDistance_)
                maxTrailingDistance_ = lastFileSize - numBytesRead;
        }

        // Try to open BAM if it wasn't opened in main loop.
        if (!bgzf) initBgzf();

        // Read any remaining data.
        while (true) {
            auto result = bgzf_read_block(bgzf.get());
            if (result != 0) {
                gziStatus_ = GziStatus::IO_ERROR;
                return;
            }
            if (bgzf->block_length == 0) break;
        }

        // Dump GZI contents to disk.
        const std::string gziFn{bamFilename_ + ".gzi"};
        ret = bgzf_index_dump(bgzf.get(), gziFn.c_str(), nullptr);
        if (ret != 0) gziStatus_ = GziStatus::GZI_ERROR;
    }

    void Write(const BamRecord& record)
    {
// TODO: add API to auto-skip this without special compile flag
#if PBBAM_AUTOVALIDATE
        Validator::Validate(record);
#endif
        // add record & its to index builder.
        //
        // NOTE: This is the record's postiion as if it were _uncompressed_. We
        //       will return with GZI data later to transform it into BAM
        //       "virtual offset".
        //
        builder_->AddRecord(record, uncompressedFilePos_);

        const auto& rawRecord = BamRecordMemory::GetRawData(record);

        // update bin
        // min_shift=14 & n_lvls=5 are BAM "magic numbers"
        rawRecord->core.bin = hts_reg2bin(rawRecord->core.pos, bam_endpos(rawRecord.get()), 14, 5);

        // write record to file
        const auto ret = sam_write1(bam_.get(), header_.get(), rawRecord.get());
        if (ret <= 0) throw IndexedBamWriterException{bamFilename_, "could not write record"};

        // update file position
        auto recordLength = [](bam1_t* b) {
            auto* c = &b->core;

            static constexpr size_t fixedLength = 36;
            const size_t qnameLength = (c->l_qname - c->l_extranul);

            // TODO: long CIGAR handling... sigh...

            size_t remainingLength = 0;
            if (c->n_cigar <= 0xffff)
                remainingLength = (b->l_data - c->l_qname);
            else {
                const size_t cigarEnd = ((uint8_t*)bam_get_cigar(b) - b->data) + (c->n_cigar * 4);
                remainingLength = 8 + (b->l_data - cigarEnd) + 4 + (4 * c->n_cigar);
            }

            return fixedLength + qnameLength + remainingLength;
        };
        uncompressedFilePos_ += recordLength(rawRecord.get());

        // Need to handle any errors from the gzi thread, since it's not set
        // up to throw without terminating the program
        auto gstatus = gziStatus_.load();
        if (gstatus != GziStatus::GOOD) {
            if (gziStatus_.load() == GziStatus::IO_ERROR)
                throw IndexedBamWriterException{bamFilename_,
                                                "error in gzi thread reading from BAM file"};
            if (gziStatus_.load() == GziStatus::TRAIL_ERROR)
                throw IndexedBamWriterException{
                    bamFilename_, "gzi reader thread failed to properly trail when reading"};
            if (gziStatus_.load() == GziStatus::GZI_ERROR)
                throw IndexedBamWriterException{bamFilename_,
                                                "could not dump GZI contents for indexing"};
            if (gziStatus_.load() == GziStatus::MISC_ERROR)
                throw IndexedBamWriterException{bamFilename_, "error computing index file"};
            gziStatus_.store(GziStatus::DEAD);
        }
    }

private:
    std::string bamFilename_;

    std::shared_ptr<bam_hdr_t> header_;
    std::unique_ptr<samFile, HtslibFileDeleter> bam_;
    std::unique_ptr<PbiBuilder2> builder_;

    // used as a type of error return code for the gziThread, so
    // that errors are delayed until at least the bam file is
    // safely written to disk
    enum class GziStatus
    {
        GOOD,
        IO_ERROR,
        TRAIL_ERROR,
        GZI_ERROR,
        MISC_ERROR,
        // There was an error, but we've bubbled up the
        // information already
        DEAD
    };
    std::atomic<GziStatus> gziStatus_{GziStatus::GOOD};
    std::thread gziThread_;

    bool isOpen_ = false;

    std::atomic<bool> done_{false};
    std::atomic<size_t> maxTrailingDistance_{0};

    int64_t uncompressedFilePos_ = 0;
};

#endif  // HTS_VERSION

IndexedBamWriter::IndexedBamWriter(const std::string& outputFilename, const BamHeader& header,
                                   const BamWriter::CompressionLevel bamCompressionLevel,
                                   const size_t numBamThreads,
                                   const PbiBuilder::CompressionLevel pbiCompressionLevel,
                                   const size_t numPbiThreads, const size_t numGziThreads,
                                   const size_t tempFileBufferSize)
    : IRecordWriter(), d_{nullptr}
{
    if (tempFileBufferSize % 8 != 0) {
        throw std::runtime_error{
            "[pbbam] indexed BAM writer ERROR: invalid buffer size for PBI builder (" +
            std::to_string(tempFileBufferSize) + "). Must be a multiple of 8."};
    }

#if PBBAM_AUTOVALIDATE
    Validator::Validate(header);
#endif
    d_ = std::make_unique<IndexedBamWriterPrivate2>(
        outputFilename, BamHeaderMemory::MakeRawHeader(header), bamCompressionLevel, numBamThreads,
        pbiCompressionLevel, numPbiThreads, numGziThreads, tempFileBufferSize);
}

IndexedBamWriter::IndexedBamWriter(IndexedBamWriter&&) noexcept = default;

IndexedBamWriter& IndexedBamWriter::operator=(IndexedBamWriter&&) noexcept = default;

IndexedBamWriter::~IndexedBamWriter() = default;

void IndexedBamWriter::Write(const BamRecord& record) { d_->Write(record); }

void IndexedBamWriter::Write(const BamRecordImpl& record) { d_->Write(BamRecord{record}); }

}  // namespace BAM
}  // namespace PacBio

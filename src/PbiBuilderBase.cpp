#include "PbbamInternalConfig.h"

#include "PbiBuilderBase.h"

#include <cassert>

namespace PacBio {
namespace BAM {

// ----------------
// PbiBuilderBase
// ----------------

PbiBuilderBase::PbiBuilderBase(const std::string& pbiFilename,
                               const PbiBuilder::CompressionLevel compressionLevel,
                               std::size_t numThreads, std::size_t bufferSize)
    : pbiFilename_{pbiFilename}
    , tempFilename_{pbiFilename + ".build"}
    , tempFile_{std::fopen(tempFilename_.c_str(), "w+b")}
    , compressionLevel_{compressionLevel}
    , numThreads_{numThreads}
    , rgIdField_{bufferSize}
    , qStartField_{bufferSize}
    , qEndField_{bufferSize}
    , holeNumField_{bufferSize}
    , readQualField_{bufferSize}
    , ctxtField_{bufferSize}
    , fileOffsetField_{bufferSize}
    , tIdField_{bufferSize}
    , tStartField_{bufferSize}
    , tEndField_{bufferSize}
    , aStartField_{bufferSize}
    , aEndField_{bufferSize}
    , revStrandField_{bufferSize}
    , nMField_{bufferSize}
    , nMMField_{bufferSize}
    , mapQualField_{bufferSize}
    , nInsOpsField_{bufferSize}
    , nDelOpsField_{bufferSize}
    , bcForwardField_{bufferSize}
    , bcReverseField_{bufferSize}
    , bcQualField_{bufferSize}
{
    if (!tempFile_) {
        throw PbiBuilderException{tempFilename_, "could not open temp file"};
    }
}

PbiBuilderBase::~PbiBuilderBase() noexcept
{
    if (!isClosed_) {
        try {
            Close();
        } catch (...) {
            // swallow any exceptions & remain no-throw from dtor
        }
    }
}

void PbiBuilderBase::AddBarcodeData(const BamRecord& b)
{
    // initialize w/ 'missing' value
    std::int16_t bcForward = -1;
    std::int16_t bcReverse = -1;
    std::int8_t bcQuality = -1;

    // check for any barcode data (both required)
    if (b.HasBarcodes() && b.HasBarcodeQuality()) {
        // fetch data from record
        std::tie(bcForward, bcReverse) = b.Barcodes();
        bcQuality = static_cast<std::int8_t>(b.BarcodeQuality());

        // double-check & reset to 'missing' value if any less than zero
        if (bcForward < 0 && bcReverse < 0 && bcQuality < 0) {
            bcForward = -1;
            bcReverse = -1;
            bcQuality = -1;
        } else {
            hasBarcodeData_ = true;
        }
    }

    // store
    bcForwardField_.Add(bcForward);
    bcReverseField_.Add(bcReverse);
    bcQualField_.Add(bcQuality);
}

void PbiBuilderBase::AddBasicData(const BamRecord& b, std::int64_t uOffset)
{
    // read group ID
    const auto rgId = [&b]() -> std::int32_t {
        auto rgIdString = b.ReadGroupBaseId();
        if (rgIdString.empty()) {
            rgIdString = MakeReadGroupId(b.MovieName(), ToString(b.Type()));
        }
        return std::stoul(rgIdString, nullptr, 16);
    }();

    // query start/end
    const auto isCcsOrTranscript = (IsCcsOrTranscript(b.Type()));
    const std::int32_t qStart = (isCcsOrTranscript ? 0 : b.QueryStart());
    const std::int32_t qEnd = (isCcsOrTranscript ? b.Impl().SequenceLength() : b.QueryEnd());

    // add'l data
    const std::int32_t holeNum = (b.HasHoleNumber() ? b.HoleNumber() : 0);
    const float readAccuracy =
        (b.HasReadAccuracy() ? boost::numeric_cast<float>(b.ReadAccuracy()) : 0.0F);
    const std::uint8_t ctxt =
        (b.HasLocalContextFlags() ? b.LocalContextFlags()
                                  : Data::LocalContextFlags::NO_LOCAL_CONTEXT);

    // store
    rgIdField_.Add(rgId);
    qStartField_.Add(qStart);
    qEndField_.Add(qEnd);
    holeNumField_.Add(holeNum);
    ctxtField_.Add(ctxt);
    readQualField_.Add(readAccuracy);
    fileOffsetField_.Add(uOffset);
}

void PbiBuilderBase::AddMappedData(const BamRecord& b)
{
    // alignment position
    const auto tId = b.ReferenceId();
    const auto tStart = static_cast<std::uint32_t>(b.ReferenceStart());
    const auto tEnd = static_cast<std::uint32_t>(b.ReferenceEnd());
    const auto aStart = static_cast<std::uint32_t>(b.AlignedStart());
    const auto aEnd = static_cast<std::uint32_t>(b.AlignedEnd());
    const auto isReverseStrand = [&b]() -> std::uint8_t {
        return (b.AlignedStrand() == Data::Strand::REVERSE ? 1 : 0);
    }();

    // alignment quality
    const auto matchData = b.NumMatchesAndMismatches();
    const auto nM = static_cast<std::uint32_t>(matchData.first);
    const auto nMM = static_cast<std::uint32_t>(matchData.second);
    const auto mapQuality = b.MapQuality();

    // indel operations
    const auto indelOps = b.NumInsertionAndDeletionOperations();
    const auto nInsOps = indelOps.first;
    const auto nDelOps = indelOps.second;

    if (tId >= 0) {
        hasMappedData_ = true;
    }

    // store
    tIdField_.Add(tId);
    tStartField_.Add(tStart);
    tEndField_.Add(tEnd);
    aStartField_.Add(aStart);
    aEndField_.Add(aEnd);
    revStrandField_.Add(isReverseStrand);
    nMField_.Add(nM);
    nMMField_.Add(nMM);
    mapQualField_.Add(mapQuality);
    nInsOpsField_.Add(nInsOps);
    nDelOpsField_.Add(nDelOps);
}

void PbiBuilderBase::AddRecord(const BamRecord& b, std::int64_t uOffset)
{
    // ensure updated data (necessary?)
    BAM::BamRecordMemory::UpdateRecordTags(b);
    b.ResetCachedPositions();

    // store record data & maybe flush to temp file
    AddBasicData(b, uOffset);
    AddMappedData(b);
    AddBarcodeData(b);
    AddReferenceData(b, currentRow_);
    FlushBuffers(FlushMode::NO_FORCE);

    ++currentRow_;
}

void PbiBuilderBase::AddReferenceData(const BamRecord& b, std::uint32_t currentRow)
{
    // only add if coordinate-sorted hint is set
    // update with info from refDataBuilder
    if (refDataBuilder_) {
        const auto sorted = refDataBuilder_->AddRecord(b, currentRow);
        if (!sorted) {
            refDataBuilder_.reset();
        }
    }
}

void PbiBuilderBase::Close()
{
    if (isClosed_) {
        return;
    }

    FlushBuffers(FlushMode::FORCE);

    OpenPbiFile();
    WritePbiHeader();
    WriteFromTempFile();

    remove(tempFilename_.c_str());
    isClosed_ = true;
}

void PbiBuilderBase::FlushBuffers(FlushMode mode)
{
    const auto force = (mode == FlushMode::FORCE);

    MaybeFlushBuffer(rgIdField_, force);
    MaybeFlushBuffer(qStartField_, force);
    MaybeFlushBuffer(qEndField_, force);
    MaybeFlushBuffer(holeNumField_, force);
    MaybeFlushBuffer(readQualField_, force);
    MaybeFlushBuffer(ctxtField_, force);
    MaybeFlushBuffer(fileOffsetField_, force);

    MaybeFlushBuffer(tIdField_, force);
    MaybeFlushBuffer(tStartField_, force);
    MaybeFlushBuffer(tEndField_, force);
    MaybeFlushBuffer(aStartField_, force);
    MaybeFlushBuffer(aEndField_, force);
    MaybeFlushBuffer(revStrandField_, force);
    MaybeFlushBuffer(nMField_, force);
    MaybeFlushBuffer(nMMField_, force);
    MaybeFlushBuffer(mapQualField_, force);
    MaybeFlushBuffer(nInsOpsField_, force);
    MaybeFlushBuffer(nDelOpsField_, force);

    MaybeFlushBuffer(bcForwardField_, force);
    MaybeFlushBuffer(bcReverseField_, force);
    MaybeFlushBuffer(bcQualField_, force);
}

void PbiBuilderBase::OpenPbiFile()
{
    // open file handle
    const auto mode = std::string("wb") + std::to_string(static_cast<int>(compressionLevel_));
    pbiFile_.reset(bgzf_open(pbiFilename_.c_str(), mode.c_str()));
    if (pbiFile_ == nullptr) {
        std::ostringstream msg;
        msg << "[pbbam] PBI index builder ERROR: could not open file for writing:\n"
            << "  file: " << pbiFilename_ << '\n';
        throw std::runtime_error{msg.str()};
    }
    // if no explicit thread count given, attempt built-in check
    std::size_t actualNumThreads = numThreads_;
    if (actualNumThreads == 0) {
        actualNumThreads = std::thread::hardware_concurrency();

        // if still unknown, default to single-threaded
        if (actualNumThreads == 0) {
            actualNumThreads = 1;
        }
    }

    // if multithreading requested, enable it
    if (actualNumThreads > 1) {
        bgzf_mt(pbiFile_.get(), actualNumThreads, 256);
    }
}

void PbiBuilderBase::WriteFromTempFile()
{
    // load from temp file, in PBI format order, and write to index

    WriteField(rgIdField_);
    WriteField(qStartField_);
    WriteField(qEndField_);
    WriteField(holeNumField_);
    WriteField(readQualField_);
    WriteField(ctxtField_);

    this->WriteVirtualOffsets();

    if (hasMappedData_) {
        WriteField(tIdField_);
        WriteField(tStartField_);
        WriteField(tEndField_);
        WriteField(aStartField_);
        WriteField(aEndField_);
        WriteField(revStrandField_);
        WriteField(nMField_);
        WriteField(nMMField_);
        WriteField(mapQualField_);
        WriteField(nInsOpsField_);
        WriteField(nDelOpsField_);
    }

    if (refDataBuilder_) {
        WriteReferenceData();
    }

    if (hasBarcodeData_) {
        WriteField(bcForwardField_);
        WriteField(bcReverseField_);
        WriteField(bcQualField_);
    }
}

void PbiBuilderBase::WritePbiHeader()
{
    BGZF* bgzf = pbiFile_.get();

    // 'magic' string
    static constexpr std::array<char, 4> MAGIC{{'P', 'B', 'I', '\1'}};
    bgzf_write_safe(bgzf, MAGIC.data(), 4);

    PbiFile::Sections sections = PbiFile::BASIC;
    if (hasMappedData_) {
        sections |= PbiFile::MAPPED;
    }
    if (hasBarcodeData_) {
        sections |= PbiFile::BARCODE;
    }
    if (refDataBuilder_) {
        sections |= PbiFile::REFERENCE;
    }

    // version, pbi_flags, & n_reads
    auto version = static_cast<std::uint32_t>(PbiFile::CurrentVersion);
    std::uint16_t pbi_flags = sections;
    auto numReads = currentRow_;
    if (bgzf->is_be) {
        version = ed_swap_4(version);
        pbi_flags = ed_swap_2(pbi_flags);
        numReads = ed_swap_4(numReads);
    }
    bgzf_write_safe(bgzf, &version, 4);
    bgzf_write_safe(bgzf, &pbi_flags, 2);
    bgzf_write_safe(bgzf, &numReads, 4);

    // reserved space
    char reserved[18];
    memset(reserved, 0, 18);
    bgzf_write_safe(bgzf, reserved, 18);
}

void PbiBuilderBase::WriteReferenceData()
{
    assert(refDataBuilder_);
    refDataBuilder_->WriteData(pbiFile_.get());
}

// -------------------------
// PbiReferenceDataBuilder
// -------------------------

PbiReferenceDataBuilder::PbiReferenceDataBuilder(std::size_t numReferenceSequences)
{
    // initialize with number of references we expect to see
    //
    // we can add more later, but want to ensure known references have an entry
    // even if no records are observed mapping to it
    //
    for (std::size_t i = 0; i < numReferenceSequences; ++i) {
        rawReferenceEntries_[i] = PbiReferenceEntry(i);
    }

    // also create an "unmapped" entry
    rawReferenceEntries_[PbiReferenceEntry::UNMAPPED_ID] = PbiReferenceEntry{};
}

bool PbiReferenceDataBuilder::AddRecord(const BamRecord& record, std::int32_t rowNumber)
{
    // fetch ref ID & pos for record
    const std::int32_t tId = record.ReferenceId();
    const std::int32_t pos = record.ReferenceStart();

    // sanity checks to protect against non-coordinate-sorted BAMs
    if (lastRefId_ != tId || (lastRefId_ >= 0 && tId < 0)) {
        if (tId >= 0) {

            // if we've already seen unmapped reads, but our current tId is valid
            //
            // error: unmapped reads should all be at the end (can stop checking refs)
            //
            PbiReferenceEntry& unmappedEntry =
                rawReferenceEntries_.at(PbiReferenceEntry::UNMAPPED_ID);
            if (unmappedEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW) {
                return false;
            }

            // if we've already seen data for this new tId
            // (remember we're coming from another tId)
            //
            // error: refs are out of order (can stop checking refs)
            //
            PbiReferenceEntry& currentEntry =
                rawReferenceEntries_.at(static_cast<std::uint32_t>(tId));
            if (currentEntry.beginRow_ != PbiReferenceEntry::UNSET_ROW) {
                return false;
            }
        }
        lastRefId_ = tId;
    } else if (tId >= 0 && lastPos_ > pos) {
        return false;  // error: positions out of order
    }

    // update row numbers
    PbiReferenceEntry& entry = rawReferenceEntries_.at(static_cast<std::uint32_t>(tId));
    if (entry.beginRow_ == PbiReferenceEntry::UNSET_ROW) {
        entry.beginRow_ = rowNumber;
    }
    entry.endRow_ = rowNumber + 1;

    // update pos (for sorting check next go-round)
    lastPos_ = pos;
    return true;
}

PbiRawReferenceData PbiReferenceDataBuilder::Result() const
{
    // PbiReferenceEntries will be sorted thanks to std::map
    // tId will be at end since we're sorting on the uint cast of -1
    PbiRawReferenceData result;
    result.entries_.reserve(rawReferenceEntries_.size());
    for (const auto& entry : rawReferenceEntries_) {
        result.entries_.push_back(entry.second);
    }
    return result;
}

void PbiReferenceDataBuilder::WriteData(BGZF* bgzf)
{
    const auto refData = Result();

    // num_refs
    std::uint32_t numRefs = refData.entries_.size();
    if (bgzf->is_be) {
        numRefs = ed_swap_4(numRefs);
    }
    bgzf_write_safe(bgzf, &numRefs, 4);

    // reference entries
    numRefs = refData.entries_.size();  // need to reset after maybe endian-swapping
    for (std::size_t i = 0; i < numRefs; ++i) {
        auto& entry = refData.entries_[i];
        auto tId = entry.tId_;
        auto beginRow = entry.beginRow_;
        auto endRow = entry.endRow_;
        if (bgzf->is_be) {
            tId = ed_swap_4(tId);
            beginRow = ed_swap_4(beginRow);
            endRow = ed_swap_4(endRow);
        }
        bgzf_write_safe(bgzf, &tId, 4);
        bgzf_write_safe(bgzf, &beginRow, 4);
        bgzf_write_safe(bgzf, &endRow, 4);
    }
}

}  // namespace BAM
}  // namespace PacBio

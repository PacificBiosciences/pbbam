#include "PbbamInternalConfig.h"

#include "ZmwChunkedFastxTextReader.h"

#include "ErrnoReason.h"

#include <memory>
#include <sstream>
#include <stdexcept>

#include <cassert>
#include <cstdio>

#include <unistd.h>

namespace PacBio {
namespace BAM {

ZmwChunkedFastxTextReader::ZmwChunkedFastxTextReader(std::string filename,
                                                     const std::size_t numChunks)
    : ZmwChunkedFastxReaderImpl{std::move(filename), numChunks}
    , file_{std::fopen(fastxFilename_.c_str(), "r")}
    , seq_{kseq_init(file_.get())}
{
    // check file handle
    if (file_ == nullptr) {
        std::ostringstream msg;
        msg << "[pbbam] chunked FASTX reader ERROR: could not open file:\n"
            << "  file: " << fastxFilename_ << '\n';
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }

    // check kseq sequence handle
    assert(seq_ != nullptr);
}

int ZmwChunkedFastxTextReader::FetchRecord(bool skipName)
{
    // NOTE: kseq_read assumes it is at the beginning of "next" sequence's name.
    //       However, here the file handle may already point to the first base after
    //       seeking using FAI. So we optionally load the name.

    int c;
    kseq_t* seq = seq_.get();
    kstream_t* ks = seq->f;
    seq_->comment.l = seq_->seq.l = seq_->qual.l = 0; /* reset all members */

    if (!skipName) {

        if (seq->last_char == 0) { /* then jump to the next header line */
            while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@') {
                ;
            }
            if (c == -1) {
                return -1; /* end of file */
            }
            seq->last_char = c;
        } /* else: the first header char has been read in the previous call */

        if (ks_getuntil(ks, 0, &seq->name, &c) < 0) {
            return -1; /* normal exit: EOF */
        }
        if (c != '\n') {
            ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */
        }
    }

    if (seq_->seq.s == 0) { /* we can do this in the loop below, but that is slower */
        seq_->seq.m = 256;
        seq_->seq.s = static_cast<char*>(std::malloc(seq_->seq.m));
    }
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        if (c == '\n') {
            continue; /* skip empty lines */
        }
        seq_->seq.s[seq_->seq.l++] = c; /* this is safe: we always have enough space for 1 char */
        ks_getuntil2(ks, KS_SEP_LINE, &seq_->seq, 0, 1); /* read the rest of the line */
    }

    if (c == '>' || c == '@') {
        seq_->last_char = c; /* the first header char has been read */
    }
    if (seq_->seq.l + 1 >=
        seq_->seq.m) { /* seq_->seq.s[seq_->seq.l] below may be out of boundary */
        seq_->seq.m = seq_->seq.l + 2;
        kroundup32(seq_->seq.m); /* rounded to the next closest 2^k */
        seq_->seq.s = static_cast<char*>(std::realloc(seq_->seq.s, seq_->seq.m));
    }
    seq_->seq.s[seq_->seq.l] = 0; /* null terminated string */

    if (c != '+') {
        return seq_->seq.l; /* FASTA */
    }
    if (seq_->qual.m < seq_->seq.m) { /* allocate memory for qual in case insufficient */
        seq_->qual.m = seq_->seq.m;
        seq_->qual.s = static_cast<char*>(std::realloc(seq_->qual.s, seq_->qual.m));
    }

    while ((c = ks_getc(ks)) != -1 && c != '\n') {
        ; /* skip the rest of '+' line */
    }
    if (c == -1) {
        return -2; /* error: no quality string */
    }
    while (ks_getuntil2(ks, KS_SEP_LINE, &seq_->qual, 0, 1) >= 0 && seq_->qual.l < seq_->seq.l) {
        ;
    }

    seq_->last_char = 0; /* we have not come to the next header line */

    if (seq_->seq.l != seq_->qual.l) {
        return -2; /* error: qual string is of a different length */
    }
    return seq_->seq.l;
}

int ZmwChunkedFastxTextReader::ReadFromFile(std::FILE* fp, void* data, std::size_t length)
{
    return static_cast<int>(std::fread(data, sizeof(std::uint8_t), length, fp));
}

FastaSequence ZmwChunkedFastxTextReader::ReadNextFasta(bool skipName)
{
    // read sequence
    const auto result = FetchRecord(skipName);
    if (result < 0) {
        std::ostringstream msg;
        msg << "[pbbam] chunked FASTX reader ERROR: could not read FASTA record:\n"
            << "  file: " << fastxFilename_ << '\n'
            << "  reason: likely truncated quality string\n";
        throw std::runtime_error{msg.str()};
    }

    // return FASTQ
    std::string name = (skipName ? "" : std::string{seq_->name.s, seq_->name.l});
    std::string bases{seq_->seq.s, seq_->seq.l};
    return FastaSequence{std::move(name), std::move(bases)};
}

FastqSequence ZmwChunkedFastxTextReader::ReadNextFastq(bool skipName)
{
    // read sequence
    const auto result = FetchRecord(skipName);
    if (result < 0) {
        std::ostringstream msg;
        msg << "[pbbam] chunked FASTX reader ERROR: could not read FASTQ record:\n"
            << "  file: " << fastxFilename_ << '\n'
            << "  reason: likely truncated quality string\n";
        throw std::runtime_error{msg.str()};
    }

    // return FASTQ
    std::string name = (skipName ? "" : std::string{seq_->name.s, seq_->name.l});
    std::string bases{seq_->seq.s, seq_->seq.l};
    Data::QualityValues quals{std::string{seq_->qual.s, seq_->qual.l}};
    return FastqSequence{std::move(name), std::move(bases), std::move(quals)};
}

void ZmwChunkedFastxTextReader::Seek(std::uint64_t pos)
{
    // seek to sequence 'id' & reset kseq handle
    auto result = std::fseek(file_.get(), pos, SEEK_SET);
    if (result != 0) {
        std::ostringstream msg;
        msg << "[pbbam] chunked FASTX reader ERROR: could not seek to requested pos: " << pos
            << '\n'
            << "  in file: " << fastxFilename_;
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }
    ks_rewind(seq_->f);
}

}  // namespace BAM
}  // namespace PacBio

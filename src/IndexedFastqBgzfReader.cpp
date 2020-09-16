#include "PbbamInternalConfig.h"

#include "IndexedFastqBgzfReader.h"

#include <cassert>

#include <algorithm>
#include <sstream>
#include <stdexcept>

#include "ErrnoReason.h"

namespace PacBio {
namespace BAM {

IndexedFastqBgzfReader::IndexedFastqBgzfReader(std::string filename)
    : IndexedFastqReaderImpl{std::move(filename)}
    , file_{bgzf_open(fastqFilename_.c_str(), "r")}
    , seq_{kseq_init(file_.get())}
{
    // check BGZF file handle
    if (file_ == nullptr) {
        std::ostringstream msg;
        msg << "[pbbam] FASTQ reader ERROR: could not open file:\n"
            << "  FASTQ file: " << fastqFilename_ << '\n';
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }

    // check kseq sequence handle
    assert(seq_ != nullptr);

    // load BGZF index data (*.gzi)
    const auto result = bgzf_index_load(file_.get(), fastqFilename_.c_str(), ".gzi");
    if (result != 0) {
        std::ostringstream msg;
        msg << "[pbbam] FASTQ reader ERROR: could not load *.gzi index data:\n"
            << "  FASTQ file: " << fastqFilename_ << '\n'
            << "  index file: " << fastqFilename_ << ".gzi";
        MaybePrintErrnoReason(msg);
        throw std::runtime_error{msg.str()};
    }
}

int IndexedFastqBgzfReader::FetchRecord()
{
    // NOTE: kseq_read assumes it is at the beginning of "next" sequence's name.
    //       However, here the file handle already points to the first base after
    //       seeking using FAI. So this is kseq_read without the name/comment scan.

    int c;
    kseq_t* seq = seq_.get();
    kstream_t* ks = seq->f;
    seq_->comment.l = seq_->seq.l = seq_->qual.l = 0; /* reset all members */
    if (seq_->seq.s == 0) { /* we can do this in the loop below, but that is slower */
        seq_->seq.m = 256;
        seq_->seq.s = (char*)malloc(seq_->seq.m);
    }
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        if (c == '\n') continue;        /* skip empty lines */
        seq_->seq.s[seq_->seq.l++] = c; /* this is safe: we always have enough space for 1 char */
        ks_getuntil2(ks, KS_SEP_LINE, &seq_->seq, 0, 1); /* read the rest of the line */
    }

    if (c == '>' || c == '@') seq_->last_char = c; /* the first header char has been read */
    if (seq_->seq.l + 1 >=
        seq_->seq.m) { /* seq_->seq.s[seq_->seq.l] below may be out of boundary */
        seq_->seq.m = seq_->seq.l + 2;
        kroundup32(seq_->seq.m); /* rounded to the next closest 2^k */
        seq_->seq.s = (char*)realloc(seq_->seq.s, seq_->seq.m);
    }
    seq_->seq.s[seq_->seq.l] = 0;     /* null terminated string */
    if (c != '+') return seq_->seq.l; /* FASTA */
    if (seq_->qual.m < seq_->seq.m) { /* allocate memory for qual in case insufficient */
        seq_->qual.m = seq_->seq.m;
        seq_->qual.s = (char*)realloc(seq_->qual.s, seq_->qual.m);
    }

    while ((c = ks_getc(ks)) != -1 && c != '\n')
        ;                   /* skip the rest of '+' line */
    if (c == -1) return -2; /* error: no quality string */
    while (ks_getuntil2(ks, KS_SEP_LINE, &seq_->qual, 0, 1) >= 0 && seq_->qual.l < seq_->seq.l)
        ;
    seq_->last_char = 0; /* we have not come to the next header line */

    if (seq_->seq.l != seq_->qual.l) return -2; /* error: qual string is of a different length */
    return seq_->seq.l;
}

std::pair<std::string, Data::QualityValues> IndexedFastqBgzfReader::Subsequence(
    const std::string& id, Data::Position start, Data::Position end)
{
    // check requested region is valid
    const auto& entry = index_.Entry(id);
    const int64_t available = static_cast<int64_t>(entry.Length) - start;
    const int64_t requested = end - start;
    const int64_t length = std::min(available, requested);
    if ((start < 0) || (end < 0) || (length < 0)) {
        std::ostringstream msg;
        msg << "[pbbam] FASTQ reader ERROR: invalid subsequence region requested:\n"
            << "  FASTQ file: " << fastqFilename_ << '\n'
            << "  requested region: " << id << ':' << start << '-' << end << '\n'
            << "  sequence length:  " << entry.Length << '\n';
        throw std::runtime_error{msg.str()};
    }

    // quick out if nothing needed
    if (length == 0) return {};

    // seek to sequence 'id' & reset kseq handle
    auto result = bgzf_useek(file_.get(), entry.SeqOffset, SEEK_SET);
    if (result != 0) {
        std::ostringstream msg;
        msg << "[pbbam] FASTQ reader ERROR: could not seek to requested region:\n"
            << "  FASTQ file: " << fastqFilename_ << '\n'
            << "  requested region: " << id << ':' << start << '-' << end << '\n';
        throw std::runtime_error{msg.str()};
    }
    ks_rewind(seq_->f);

    // read (entire) sequence
    result = FetchRecord();
    if (result < 0) {
        std::ostringstream msg;
        msg << "[pbbam] FASTQ reader ERROR: could not read FASTQ record:\n"
            << "  FASTQ file: " << fastqFilename_ << '\n'
            << "  requested region: " << id << ':' << start << '-' << end << '\n'
            << "  reason: likely truncated quality string\n";
        throw std::runtime_error{msg.str()};
    }

    // trim to region bounds and return
    const std::string seq{seq_->seq.s, seq_->seq.l};
    const std::string quals{seq_->qual.s, seq_->qual.l};
    return std::make_pair(seq.substr(start, length),
                          Data::QualityValues::FromFastq(quals.substr(start, length)));
}

}  // namespace BAM
}  // namespace PacBio

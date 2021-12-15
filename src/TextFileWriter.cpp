#include "PbbamInternalConfig.h"

#include <pbbam/TextFileWriter.h>

#include <pbbam/Deleters.h>
#include "ErrnoReason.h"
#include "FileProducer.h"

#include <boost/algorithm/string.hpp>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace PacBio {
namespace BAM {

class TextFileWriter::TextFileWriterPrivate : public FileProducer
{
public:
    TextFileWriterPrivate(const std::string& filename) : FileProducer{filename}
    {
        isZipped_ = boost::algorithm::iends_with(filename, ".gz");

        if (isZipped_) {
            // open for gzipped text
            bgzf_.reset(bgzf_open(TempFilename().c_str(), "wg"));
            if (bgzf_.get() == nullptr) {
                std::ostringstream msg;
                msg << "[pbbam] text file writer ERROR: could not open zipped file:\n"
                    << "  file: " << filename;
                MaybePrintErrnoReason(msg);
                throw std::runtime_error{msg.str()};
            }
        } else {
            // open for plain text
            out_.open(TempFilename());
            if (!out_) {
                std::ostringstream msg;
                msg << "[pbbam] text file writer ERROR: could not open plain text file:\n"
                    << "  file: " << filename;
                MaybePrintErrnoReason(msg);
                throw std::runtime_error{msg.str()};
            }
        }
    }

    void Write(const std::string& line)
    {
        if (isZipped_) {
            const size_t length = line.size();
            ssize_t written = bgzf_write(bgzf_.get(), line.c_str(), length);
            written += bgzf_write(bgzf_.get(), "\n", 1);
            if (written != static_cast<ssize_t>(length + 1)) {
                std::ostringstream msg;
                msg << "[pbbam] text file writer ERROR: could not write to file:\n"
                    << "  file: " << TempFilename();
                MaybePrintErrnoReason(msg);
                throw std::runtime_error{msg.str()};
            }
        } else {
            out_ << line << '\n';
        }
    }

    bool isZipped_ = false;
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf_;
    std::ofstream out_;
};

TextFileWriter::TextFileWriter(const std::string& filename)
    : d_{std::make_unique<TextFileWriterPrivate>(filename)}
{
}

TextFileWriter::TextFileWriter(TextFileWriter&&) noexcept = default;

TextFileWriter& TextFileWriter::operator=(TextFileWriter&&) noexcept = default;

TextFileWriter::~TextFileWriter() = default;

void TextFileWriter::Write(const std::string& line) { d_->Write(line); }

}  // namespace BAM
}  // namespace PacBio


// File Description
/// \file TextFileWriter.cpp
/// \brief Implements the TextFileWriter class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/TextFileWriter.h"

#include <cassert>

#include <fstream>
#include <iostream>
#include <type_traits>

#include <boost/algorithm/string.hpp>

#include "FileProducer.h"
#include "MemoryUtils.h"

namespace PacBio {
namespace BAM {

static_assert(!std::is_copy_constructible<TextFileWriter>::value,
              "TextFileWriter(const TextFileWriter&) is not = delete");
static_assert(!std::is_copy_assignable<TextFileWriter>::value,
              "TextFileWriter& operator=(const TextFileWriter&) is not = delete");

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
                throw std::runtime_error("TextFileWriter - could not open file: " + filename +
                                         " for writing");
            }
        } else {
            // open for plain text
            out_.open(TempFilename());
            if (!out_) {
                throw std::runtime_error("TextFileWriter - could not open file: " + filename +
                                         " for writing");
            }
        }
    }

    void Write(const std::string& line)
    {
        if (isZipped_) {
            const size_t length = line.size();
            ssize_t written = bgzf_write(bgzf_.get(), line.c_str(), length);
            written += bgzf_write(bgzf_.get(), "\n", 1);
            if (written != static_cast<ssize_t>(length + 1))
                throw std::runtime_error("TextFileWriter - error writing to file: " +
                                         TargetFilename());
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
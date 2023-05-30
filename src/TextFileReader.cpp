#include "PbbamInternalConfig.h"

#include <pbbam/TextFileReader.h>

#include <pbbam/Deleters.h>
#include "ErrnoReason.h"

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace PacBio {
namespace BAM {

class TextFileReader::TextFileReaderPrivate
{
public:
    TextFileReaderPrivate(std::string filename) : filename_{std::move(filename)}, k_{0, 0, nullptr}
    {
        // "ru" : read & supply plain output
        bgzf_.reset(bgzf_open(filename_.c_str(), "ru"));
        if (bgzf_.get() == nullptr) {
            std::ostringstream msg;
            msg << "[pbbam] text file reader ERROR: could not open file:\n"
                << "  file: " << filename_;
            MaybePrintErrnoReason(msg);
            throw std::runtime_error{msg.str()};
        }

        // pre-fetch first line
        GetNext();
    }

    ~TextFileReaderPrivate() { std::free(k_.s); }

    void GetNext()
    {
        line_.clear();

        // be sure we skip empty lines
        while (line_.empty()) {
            const int result = bgzf_getline(bgzf_.get(), '\n', &k_);

            // found data
            if (result > 0) {
                line_ = std::string{k_.s, k_.l};
                return;
            }

            // empty line, try again
            else if (result == 0) {
                continue;

                // EOF (no error, but will stop next TextFileReader iteration
            } else if (result == -1) {
                return;

                // else error
            } else {
                std::ostringstream msg;
                msg << "[pbbam] text file reader ERROR: could not read from file:\n"
                    << "  file: " << filename_ << '\n'
                    << "  htslib error code: " << result;
                MaybePrintErrnoReason(msg);
                throw std::runtime_error{msg.str()};
            }
        }
    }

    std::string filename_;
    std::string line_;
    std::unique_ptr<BGZF, HtslibBgzfDeleter> bgzf_;  // can handle plain text or gzipped
    kstring_t k_;
};

TextFileReader::TextFileReader(std::string filename)
    : PacBio::BAM::internal::QueryBase<std::string>{}
    , d_{std::make_unique<TextFileReaderPrivate>(std::move(filename))}
{}

TextFileReader::TextFileReader(TextFileReader&&) noexcept = default;

TextFileReader& TextFileReader::operator=(TextFileReader&&) noexcept = default;

TextFileReader::~TextFileReader() = default;

const std::string& TextFileReader::Filename() const { return d_->filename_; }

bool TextFileReader::GetNext(std::string& line)
{
    if (d_->line_.empty()) {
        return false;
    }

    line = d_->line_;
    d_->GetNext();
    return true;
}

std::vector<std::string> TextFileReader::ReadAll(const std::string& fn)
{
    std::vector<std::string> result;
    result.reserve(256);
    TextFileReader reader{fn};
    for (const auto& seq : reader) {
        result.emplace_back(seq);
    }
    return result;
}

}  // namespace BAM
}  // namespace PacBio

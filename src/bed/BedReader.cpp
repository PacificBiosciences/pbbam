#include "PbbamInternalConfig.h"

#include <pbbam/bed/BedReader.h>

#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include <pbbam/FormatUtils.h>
#include <pbbam/StringUtilities.h>
#include <pbbam/TextFileReader.h>

namespace PacBio {
namespace BED {

class BedReader::BedReaderPrivate
{
public:
    explicit BedReaderPrivate(const std::string& fn)
    {
        // validate extension
        if (!BAM::FormatUtils::IsBedFilename(fn)) {
            std::ostringstream msg;
            msg << "[pbbam] BED reader ERROR: not a recognized BED extension:\n"
                << "  filename: " << fn << '\n';
            throw std::runtime_error{msg.str()};
        };

        // open file stream
        reader_ = std::make_unique<BAM::TextFileReader>(fn);
        if (!reader_) {
            std::ostringstream msg;
            msg << "[pbbam] BED reader ERROR: could not open file:\n"
                << "  BED file: " << fn << '\n';
            throw std::runtime_error{msg.str()};
        }

        // pre-fetch first record
        GetNext();
    }

    const std::string& Filename() const { return reader_->Filename(); }

    void GetNext()
    {
        interval_ = boost::none;
        std::string line;
        if (reader_->GetNext(line)) {
            interval_ = ParseInterval(std::move(line));
        }
    }

    Data::GenomicInterval ParseInterval(std::string line)
    {
        // trim any trailing whitespace
        boost::trim_right(line);

        // split into token fields
        const auto fields = BAM::Split(line, '\t');
        if (fields.size() < 3) {
            std::ostringstream msg;
            msg << "[pbbam] BED reader ERROR: invalid BED record. Line has fewer than 3 fields:\n"
                << line << "'\n"
                << "  BED file: " << reader_->Filename() << '\n';
            throw std::runtime_error{msg.str()};
        }

        // convert fields into interval
        const Data::Position start = std::stoi(fields[1]);
        const Data::Position end = std::stoi(fields[2]);
        return {fields[0], start, end};
    }

    std::unique_ptr<BAM::TextFileReader> reader_;
    boost::optional<Data::GenomicInterval> interval_;
};

BedReader::BedReader(const std::string& fn)
    : BAM::internal::QueryBase<Data::GenomicInterval>{}, d_{std::make_unique<BedReaderPrivate>(fn)}
{
}

BedReader::BedReader(BedReader&&) noexcept = default;

BedReader& BedReader::operator=(BedReader&&) noexcept = default;

BedReader::~BedReader() = default;

const std::string& BedReader::Filename() const { return d_->Filename(); }

bool BedReader::GetNext(Data::GenomicInterval& interval)
{
    if (!d_->interval_) {
        return false;
    }

    interval = *d_->interval_;
    d_->GetNext();
    return true;
}

std::vector<Data::GenomicInterval> BedReader::ReadAll(const std::string& fn)
{
    std::vector<Data::GenomicInterval> result;
    result.reserve(256);
    BedReader reader{fn};
    for (const auto& seq : reader) {
        result.emplace_back(seq);
    }
    return result;
}

}  // namespace BED
}  // namespace PacBio

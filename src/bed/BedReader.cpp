// File Description
/// \file BedReader.cpp
/// \brief Implements the BedReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/bed/BedReader.h"

#include <cassert>

#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "pbbam/FormatUtils.h"
#include "pbbam/StringUtilities.h"
#include "pbbam/TextFileReader.h"

namespace PacBio {
namespace BAM {

static_assert(!std::is_copy_constructible<BedReader>::value,
              "BedReader(const BedReader&) is not = delete");
static_assert(!std::is_copy_assignable<BedReader>::value,
              "BedReader& operator=(const BedReader&) is not = delete");

class BedReader::BedReaderPrivate
{
public:
    explicit BedReaderPrivate(const std::string& fn)
    {
        // validate extension
        if (!FormatUtils::IsBedFilename(fn)) {
            throw std::runtime_error{"BedReader ERROR: filename '" + fn +
                                     "' is not recognized as a BED file."};
        }

        // open file stream
        reader_ = std::make_unique<TextFileReader>(fn);
        if (!reader_) {
            throw std::runtime_error("BedReader ERROR: could not open text file '" + fn +
                                     "' for reading");
        }

        // pre-fetch first record
        GetNext();
    }

    void GetNext()
    {
        interval_ = boost::none;
        std::string line;
        if (reader_->GetNext(line)) interval_ = ParseInterval(std::move(line));
    }

    GenomicInterval ParseInterval(std::string line)
    {
        // trim any trailing whitespace
        boost::trim_right(line);

        // split into token fields
        const auto fields = PacBio::BAM::Split(line, '\t');
        if (fields.size() < 3) {
            std::ostringstream msg;
            msg << "BedReader ERROR: invalid BED record. Line:\n"
                << line << '\n'
                << "has less than 3 fields.";
            throw std::runtime_error{msg.str()};
        }

        // convert fields into interval
        const Position start = std::stoi(fields[1]);
        const Position end = std::stoi(fields[2]);
        return {fields[0], start, end};
    }

    std::unique_ptr<TextFileReader> reader_;
    boost::optional<GenomicInterval> interval_;
};

BedReader::BedReader(const std::string& fn)
    : internal::QueryBase<GenomicInterval>{}, d_{std::make_unique<BedReaderPrivate>(fn)}
{
}

BedReader::BedReader(BedReader&&) noexcept = default;

BedReader& BedReader::operator=(BedReader&&) noexcept = default;

BedReader::~BedReader() = default;

bool BedReader::GetNext(GenomicInterval& interval)
{
    if (!d_->interval_) return false;

    interval = *d_->interval_;
    d_->GetNext();
    return true;
}

std::vector<GenomicInterval> BedReader::ReadAll(const std::string& fn)
{
    std::vector<GenomicInterval> result;
    result.reserve(256);
    BedReader reader{fn};
    for (const auto& seq : reader)
        result.emplace_back(seq);
    return result;
}

}  // namespace BAM
}  // namespace PacBio

// File Description
/// \file FastaReader.cpp
/// \brief Implements the FastaReader class.
//
// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/FaiIndex.h"

#include <cassert>

#include <fstream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <htslib/faidx.h>

#include "pbbam/StringUtilities.h"

namespace PacBio {
namespace BAM {

static_assert(!std::is_copy_constructible<FaiIndex>::value,
              "FaiIndex(const FaiIndex&) is not = delete");
static_assert(!std::is_copy_assignable<FaiIndex>::value,
              "FaiIndex& operator=(const FaiIndex&) is not = delete");

bool operator==(const FaiEntry& lhs, const FaiEntry& rhs) noexcept
{
    return std::tie(lhs.Length, lhs.SeqOffset, lhs.NumBases, lhs.NumBytes, lhs.QualOffset) ==
           std::tie(rhs.Length, rhs.SeqOffset, rhs.NumBases, rhs.NumBytes, rhs.QualOffset);
}

std::ostream& operator<<(std::ostream& out, const FaiEntry& entry)
{
    out << entry.Length << '\t' << entry.SeqOffset << '\t' << entry.NumBases << '\t'
        << entry.NumBytes;
    if (entry.QualOffset >= 0) out << '\t' << entry.QualOffset;
    return out;
}

class FaiIndex::FaiIndexPrivate
{
public:
    FaiIndexPrivate() = default;
    FaiIndexPrivate(const std::string& fn) { LoadFromFile(fn); }

    void Add(std::string name, FaiEntry entry)
    {
        names_.push_back(name);
        data_.emplace(std::move(name), std::move(entry));
    }

    void LoadFromFile(const std::string& fn)
    {
        std::ifstream f{fn};
        std::string line;
        std::vector<std::string> fields;
        while (std::getline(f, line)) {

            fields = Split(line, '\t');
            const auto numFields = fields.size();
            if (numFields < 5 || numFields > 6) {
                std::ostringstream msg;
                msg << "[pbbam] FAI index ERROR: malformed index line, incorrect number of "
                       "fields:\n"
                    << "  expected: 5 for FASTA, or 6 for FASTQ\n"
                    << "  observed: " << numFields << " in line:\n"
                    << line << '\n'
                    << "  file: " << fn;
                throw std::runtime_error{msg.str()};
            }

            FaiEntry entry;
            entry.Length = std::stoull(fields[1]);
            entry.SeqOffset = std::stoull(fields[2]);
            entry.NumBases = std::stoul(fields[3]);
            entry.NumBytes = std::stoul(fields[4]);
            if (numFields == 6) entry.QualOffset = std::stoll(fields[5]);

            Add(std::move(fields[0]), std::move(entry));
        }
    }

    std::vector<std::string> names_;                  // save names in input order
    std::unordered_map<std::string, FaiEntry> data_;  // map name -> data
};

FaiIndex::FaiIndex(const std::string& fn) : d_{std::make_unique<FaiIndexPrivate>(fn)} {}

FaiIndex::FaiIndex() : d_{std::make_unique<FaiIndexPrivate>()} {}

FaiIndex::FaiIndex(FaiIndex&&) noexcept = default;

FaiIndex& FaiIndex::operator=(FaiIndex&&) noexcept = default;

FaiIndex::~FaiIndex() = default;

void FaiIndex::Add(std::string name, FaiEntry entry) { d_->Add(std::move(name), std::move(entry)); }

void FaiIndex::Create(const std::string& fn)
{
    const auto ret = fai_build(fn.c_str());
    if (ret < 0) {
        throw std::runtime_error{"[pbbam] FAI index ERROR: could not create *.fai for file:\n" +
                                 fn};
    }
}

const FaiEntry& FaiIndex::Entry(const std::string& name) const
{
    const auto found = d_->data_.find(name);
    if (found == d_->data_.cend())
        throw std::runtime_error{
            "[pbbam] FAI index ERROR: could not find entry for sequence name: " + name};
    return found->second;
}

const FaiEntry& FaiIndex::Entry(const uint32_t row) const
{
    const auto& name = d_->names_.at(row);
    return Entry(name);
}

bool FaiIndex::HasEntry(const std::string& name) const
{
    const auto found = d_->data_.find(name);
    return found != d_->data_.cend();
}

const std::vector<std::string>& FaiIndex::Names() const { return d_->names_; }

void FaiIndex::Save(const std::string& fn) const
{
    std::ofstream out{fn};
    Save(out);
}

void FaiIndex::Save(std::ostream& out) const
{
    for (const auto& name : d_->names_) {
        const auto& entry = Entry(name);
        out << name << '\t' << entry << '\n';
    }
}

}  // namespace BAM
}  // namespace PacBio

#include "PbbamInternalConfig.h"

#include <pbbam/BamHeader.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <htslib/hts.h>

#include <pbbam/BamFile.h>
#include <pbbam/DataSet.h>
#include <pbbam/SamTagCodec.h>
#include <pbbam/StringUtilities.h>

#include "Version.h"

namespace PacBio {
namespace BAM {
namespace {

const std::string BamHeaderPrefixHD{"@HD"};
const std::string BamHeaderPrefixSQ{"@SQ"};
const std::string BamHeaderPrefixRG{"@RG"};
const std::string BamHeaderPrefixPG{"@PG"};
const std::string BamHeaderPrefixCO{"@CO"};

const std::string BamHeaderTokenVN{"VN"};
const std::string BamHeaderTokenSO{"SO"};
const std::string BamHeaderTokenpb{"pb"};

const std::string CurrentSamFormatVersion{"1.6"};

bool CheckSortOrder(const std::string& lhs, const std::string& rhs) { return lhs == rhs; }

bool CheckPbVersion(const std::string& lhs, const std::string& rhs)
{
    return (Version{lhs} >= Version::Minimum && Version{rhs} >= Version::Minimum);
}

bool CheckSequences(const std::string& sortOrder, const std::vector<SequenceInfo>& lhs,
                    const std::vector<SequenceInfo>& rhs)
{
    return ((sortOrder == "coordinate") ? lhs == rhs : true);
}

static void EnsureCanMerge(const BamHeader& lhs, const BamHeader& rhs)
{
    // check compatibility
    const auto sortOrderOk = CheckSortOrder(lhs.SortOrder(), rhs.SortOrder());
    const auto pbVersionOk = CheckPbVersion(lhs.PacBioBamVersion(), rhs.PacBioBamVersion());
    const auto sequencesOk = CheckSequences(lhs.SortOrder(), lhs.Sequences(), rhs.Sequences());
    if (sortOrderOk && pbVersionOk && sequencesOk) {
        return;
    }

    // if any checks failed, format error message & throw
    std::ostringstream e;
    e << "[pbbam] BAM header ERROR: could not merge headers:\n";

    if (!sortOrderOk) {
        e << "  mismatched sort orders (@HD:SO) : (" << lhs.SortOrder() << ", " << rhs.SortOrder()
          << ")\n";
    }

    if (!pbVersionOk) {
        e << "  incompatible PacBio BAM versions (@HD:pb) : (" << lhs.PacBioBamVersion() << ", "
          << rhs.PacBioBamVersion() << ")\n";
    }

    if (!sequencesOk) {
        e << "  mismatched sequence lists (@SQ entries)\n";
    }

    throw std::runtime_error{e.str()};
}

void ParseHeaderLine(const std::string& line, BamHeader& hdr)
{
    // pop off '@HD\t', then split HD lines into tokens
    const auto tokens = Split(line.substr(4), '\t');
    for (const auto& token : tokens) {
        const auto tokenTag = token.substr(0, 2);
        auto tokenValue = token.substr(3);

        // set header contents
        if (tokenTag == BamHeaderTokenVN) {
            hdr.Version(std::move(tokenValue));
        } else if (tokenTag == BamHeaderTokenSO) {
            hdr.SortOrder(std::move(tokenValue));
        } else if (tokenTag == BamHeaderTokenpb) {
            hdr.PacBioBamVersion(std::move(tokenValue));
        }
    }
}

}  // namespace

class BamHeader::BamHeaderPrivate
{
public:
    std::string version_;
    std::string pacbioBamVersion_;
    std::string sortOrder_;
    std::map<std::string, std::string> headerLineCustom_;

    std::map<std::string, ReadGroupInfo> readGroups_;  // id => read group info
    std::map<std::string, ProgramInfo> programs_;      // id => program info
    std::vector<std::string> comments_;

    // we need to preserve insertion order, use lookup for access by name
    std::vector<SequenceInfo> sequences_;
    std::map<std::string, int32_t> sequenceIdLookup_;
};

BamHeader::BamHeader() : d_{std::make_shared<BamHeaderPrivate>()} {}

BamHeader::BamHeader(const DataSet& dataset) : BamHeader{dataset.BamFilenames()} {}

BamHeader::BamHeader(const std::vector<std::string>& bamFilenames) : BamHeader{}
{
    if (bamFilenames.empty()) {
        throw std::runtime_error{"[pbbam] BAM header merging ERROR: no input filenames provided"};
    }

    std::vector<BamHeader> headers;
    for (const auto& fn : bamFilenames) {
        const BamFile bamFile{fn};
        headers.push_back(bamFile.Header());
    }

    *this = headers.at(0);
    for (size_t i = 1; i < headers.size(); ++i) {
        *this += headers.at(i);
    }
}

BamHeader::BamHeader(const std::vector<BamHeader>& headers) : BamHeader()
{
    if (headers.empty()) {
        throw std::runtime_error{"[pbbam] BAM header merging ERROR: no input headers provided"};
    }

    *this = headers.at(0);
    for (size_t i = 1; i < headers.size(); ++i) {
        *this += headers.at(i);
    }
}

BamHeader::BamHeader(const std::string& samHeaderText) : d_{std::make_shared<BamHeaderPrivate>()}
{
    std::istringstream s{samHeaderText};
    std::string line;
    std::string firstToken;
    while (std::getline(s, line)) {

        // skip if line is not long enough to contain true values
        if (line.length() < 5) {
            continue;
        }

        // determine token at beginning of line
        firstToken = line.substr(0, 3);

        if (firstToken == BamHeaderPrefixHD) {
            ParseHeaderLine(line, *this);
            if (Version().empty()) {
                Version(std::string{hts_version()});
            }
        }

        else if (firstToken == BamHeaderPrefixSQ) {
            AddSequence(SequenceInfo::FromSam(line));

        } else if (firstToken == BamHeaderPrefixRG) {
            AddReadGroup(ReadGroupInfo::FromSam(line));

        } else if (firstToken == BamHeaderPrefixPG) {
            AddProgram(ProgramInfo::FromSam(line));

        } else if (firstToken == BamHeaderPrefixCO) {
            AddComment(line.substr(4));
        }
    }
}

BamHeader& BamHeader::operator+=(const BamHeader& other)
{
    EnsureCanMerge(*this, other);

    // merge read groups
    for (const auto& rg : other.ReadGroups()) {
        if (!HasReadGroup(rg.Id())) {
            AddReadGroup(rg);
        }
    }

    // merge programs
    for (const auto& pg : other.Programs()) {
        if (!HasProgram(pg.Id())) {
            AddProgram(pg);
        }
    }

    // merge comments
    for (const auto& comment : other.Comments()) {
        AddComment(comment);
    }

    return *this;
}

BamHeader BamHeader::operator+(const BamHeader& other) const { return DeepCopy() += other; }

BamHeader& BamHeader::AddComment(std::string comment)
{
    d_->comments_.push_back(std::move(comment));
    return *this;
}

BamHeader& BamHeader::AddProgram(ProgramInfo pg)
{
    d_->programs_[pg.Id()] = std::move(pg);
    return *this;
}

BamHeader& BamHeader::AddReadGroup(ReadGroupInfo readGroup)
{
    const auto id = readGroup.Id();
    if (!HasReadGroup(id)) {
        d_->readGroups_[id] = std::move(readGroup);
    }
    return *this;
}

BamHeader& BamHeader::AddSequence(SequenceInfo sequence)
{
    const auto name = sequence.Name();
    if (!HasSequence(name)) {
        d_->sequences_.push_back(std::move(sequence));
        d_->sequenceIdLookup_[name] = d_->sequences_.size() - 1;
    }
    return *this;
}

BamHeader& BamHeader::ClearComments()
{
    d_->comments_.clear();
    return *this;
}

BamHeader& BamHeader::ClearPrograms()
{
    d_->programs_.clear();
    return *this;
}

BamHeader& BamHeader::ClearReadGroups()
{
    d_->readGroups_.clear();
    return *this;
}

BamHeader& BamHeader::ClearSequences()
{
    d_->sequenceIdLookup_.clear();
    d_->sequences_.clear();
    return *this;
}

std::vector<std::string> BamHeader::Comments() const { return d_->comments_; }

BamHeader& BamHeader::Comments(std::vector<std::string> comments)
{
    d_->comments_ = std::move(comments);
    return *this;
}

BamHeader BamHeader::DeepCopy() const
{
    BamHeader result;
    result.d_->version_ = d_->version_;
    result.d_->pacbioBamVersion_ = d_->pacbioBamVersion_;
    result.d_->sortOrder_ = d_->sortOrder_;
    result.d_->headerLineCustom_ = d_->headerLineCustom_;
    result.d_->readGroups_ = d_->readGroups_;
    result.d_->programs_ = d_->programs_;
    result.d_->comments_ = d_->comments_;
    result.d_->sequences_ = d_->sequences_;
    result.d_->sequenceIdLookup_ = d_->sequenceIdLookup_;
    return result;
}

bool BamHeader::HasProgram(const std::string& id) const
{
    return d_->programs_.find(id) != d_->programs_.cend();
}

bool BamHeader::HasReadGroup(const std::string& id) const
{
    return d_->readGroups_.find(id) != d_->readGroups_.cend();
}

bool BamHeader::HasSequence(const std::string& name) const
{
    return d_->sequenceIdLookup_.find(name) != d_->sequenceIdLookup_.cend();
}

size_t BamHeader::NumSequences() const { return d_->sequences_.size(); }

bool BamHeader::Empty() const noexcept
{
    assert(d_);
    return d_->version_.empty() && d_->pacbioBamVersion_.empty() && d_->sortOrder_.empty() &&
           d_->headerLineCustom_.empty() && d_->readGroups_.empty() && d_->programs_.empty() &&
           d_->comments_.empty() && d_->sequences_.empty() && d_->sequenceIdLookup_.empty();
}

std::string BamHeader::PacBioBamVersion() const { return d_->pacbioBamVersion_; }

BamHeader& BamHeader::PacBioBamVersion(const std::string& version)
{
    d_->pacbioBamVersion_ = version;
    const BAM::Version fileVersion{version};
    if (fileVersion < Version::Minimum) {
        throw std::runtime_error{"[pbbam] BAM header ERROR: invalid PacBio BAM version number (" +
                                 fileVersion.ToString() +
                                 ") is older than the minimum supported version (" +
                                 Version::Minimum.ToString() + ")"};
    }
    return *this;
}

ProgramInfo BamHeader::Program(const std::string& id) const
{
    const auto iter = d_->programs_.find(id);
    if (iter == d_->programs_.cend()) {
        throw std::runtime_error{"[pbbam] BAM header ERROR: program ID not found: " + id};
    }
    return iter->second;
}

std::vector<std::string> BamHeader::ProgramIds() const
{
    std::vector<std::string> result;
    result.reserve(d_->programs_.size());
    for (const auto& pg : d_->programs_) {
        result.push_back(pg.first);
    }
    return result;
}

std::vector<ProgramInfo> BamHeader::Programs() const
{
    std::vector<ProgramInfo> result;
    result.reserve(d_->programs_.size());
    for (const auto& pg : d_->programs_) {
        result.push_back(pg.second);
    }
    return result;
}

BamHeader& BamHeader::Programs(std::vector<ProgramInfo> programs)
{
    d_->programs_.clear();
    for (const auto& pg : programs) {
        d_->programs_[pg.Id()] = std::move(pg);
    }
    return *this;
}

ReadGroupInfo BamHeader::ReadGroup(const std::string& id) const
{
    const auto iter = d_->readGroups_.find(id);
    if (iter == d_->readGroups_.cend()) {
        throw std::runtime_error{"[pbbam] BAM header ERROR: read group ID not found: " + id};
    }
    return iter->second;
}

std::vector<std::string> BamHeader::ReadGroupIds() const
{
    std::vector<std::string> result;
    result.reserve(d_->readGroups_.size());
    for (const auto& rg : d_->readGroups_) {
        result.push_back(rg.first);
    }
    return result;
}

std::vector<ReadGroupInfo> BamHeader::ReadGroups() const
{
    std::vector<ReadGroupInfo> result;
    result.reserve(d_->readGroups_.size());
    for (const auto& rg : d_->readGroups_) {
        result.push_back(rg.second);
    }
    return result;
}

BamHeader& BamHeader::ReadGroups(std::vector<ReadGroupInfo> readGroups)
{
    d_->readGroups_.clear();
    for (auto&& rg : readGroups) {
        AddReadGroup(std::move(rg));
    }
    return *this;
}

SequenceInfo BamHeader::Sequence(const int32_t id) const { return d_->sequences_.at(id); }

SequenceInfo BamHeader::Sequence(const std::string& name) const
{
    // TODO: SequenceId(name) throws if not found, should we do so here as well?

    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend()) {
        return SequenceInfo();
    }
    const auto index = iter->second;
    assert(index >= 0 && static_cast<size_t>(index) < d_->sequences_.size());
    return d_->sequences_.at(index);
}

int32_t BamHeader::SequenceId(const std::string& name) const
{
    const auto iter = d_->sequenceIdLookup_.find(name);
    if (iter == d_->sequenceIdLookup_.cend()) {
        throw std::runtime_error{"[pbbam] BAM header ERROR: sequence name not found: " + name};
    }
    return iter->second;
}

std::string BamHeader::SequenceLength(const int32_t id) const { return Sequence(id).Length(); }

std::string BamHeader::SequenceName(const int32_t id) const { return Sequence(id).Name(); }

std::vector<std::string> BamHeader::SequenceNames() const
{
    std::vector<std::string> result;
    result.reserve(d_->sequences_.size());
    for (const auto& seq : d_->sequences_) {
        result.push_back(seq.Name());
    }
    return result;
}

std::vector<SequenceInfo> BamHeader::Sequences() const { return d_->sequences_; }

BamHeader& BamHeader::Sequences(std::vector<SequenceInfo> sequences)
{
    d_->sequences_.clear();
    for (auto&& seq : sequences) {
        AddSequence(std::move(seq));
    }
    return *this;
}

std::string BamHeader::SortOrder() const { return d_->sortOrder_; }

BamHeader& BamHeader::SortOrder(std::string order)
{
    d_->sortOrder_ = std::move(order);
    return *this;
}

std::string BamHeader::ToSam() const
{
    // init stream
    std::ostringstream out;

    // @HD
    const auto outputVersion = (d_->version_.empty() ? CurrentSamFormatVersion : d_->version_);
    const auto outputSortOrder = (d_->sortOrder_.empty() ? std::string{"unknown"} : d_->sortOrder_);
    const auto outputPbBamVersion =
        (d_->pacbioBamVersion_.empty() ? Version::Current.ToString() : d_->pacbioBamVersion_);

    out << BamHeaderPrefixHD << MakeSamTag(BamHeaderTokenVN, outputVersion)
        << MakeSamTag(BamHeaderTokenSO, outputSortOrder)
        << MakeSamTag(BamHeaderTokenpb, outputPbBamVersion) << '\n';

    // @SQ
    for (const auto& seq : d_->sequences_) {
        out << seq.ToSam() << '\n';
    }

    // @RG
    for (const auto& rgIter : d_->readGroups_) {
        out << rgIter.second.ToSam() << '\n';
    }

    // @PG
    for (const auto& progIter : d_->programs_) {
        out << progIter.second.ToSam() << '\n';
    }

    // @CO
    for (const auto& comment : d_->comments_) {
        out << BamHeaderPrefixCO << '\t' << comment << '\n';
    }

    // return result
    return out.str();
}

std::string BamHeader::Version() const { return d_->version_; }

BamHeader& BamHeader::Version(std::string version)
{
    d_->version_ = std::move(version);
    return *this;
}

}  // namespace BAM
}  // namespace PacBio

// File Description
/// \file BamHeader.inl
/// \brief Inline implementations for the BamHeader class.
//
// Author: Derek Barnett

#include "pbbam/BamHeader.h"

namespace PacBio {
namespace BAM {
namespace internal {

class BamHeaderPrivate
{
public:
    std::string version_;
    std::string pacbioBamVersion_;
    std::string sortOrder_;
    std::map<std::string, std::string> headerLineCustom_;

    std::map<std::string, ReadGroupInfo> readGroups_; // id => read group info
    std::map<std::string, ProgramInfo> programs_;     // id => program info
    std::vector<std::string> comments_;

    // we need to preserve insertion order, use lookup for access by name
    std::vector<SequenceInfo> sequences_;
    std::map<std::string, int32_t> sequenceIdLookup_;
};

} // namespace internal

inline BamHeader::BamHeader()
    : d_{std::make_shared<internal::BamHeaderPrivate>()}
{ }

inline BamHeader BamHeader::operator+(const BamHeader& other) const
{ return DeepCopy() += other; }

inline BamHeader& BamHeader::AddComment(std::string comment)
{ d_->comments_.push_back(std::move(comment)); return *this; }

inline BamHeader& BamHeader::AddProgram(ProgramInfo pg)
{ d_->programs_[pg.Id()] = std::move(pg); return *this; }

inline BamHeader& BamHeader::AddReadGroup(ReadGroupInfo readGroup)
{ d_->readGroups_[readGroup.Id()] = std::move(readGroup); return *this; }

inline BamHeader& BamHeader::ClearComments()
{ d_->comments_.clear(); return* this; }

inline BamHeader& BamHeader::ClearPrograms()
{ d_->programs_.clear(); return *this; }

inline BamHeader& BamHeader::ClearReadGroups()
{ d_->readGroups_.clear(); return *this; }

inline std::vector<std::string> BamHeader::Comments() const
{ return d_->comments_; }

inline BamHeader& BamHeader::Comments(std::vector<std::string> comments)
{ d_->comments_ = std::move(comments); return *this; }

inline bool BamHeader::HasProgram(const std::string& id) const
{ return d_->programs_.find(id) != d_->programs_.cend(); }

inline bool BamHeader::HasReadGroup(const std::string& id) const
{ return d_->readGroups_.find(id) != d_->readGroups_.cend(); }

inline bool BamHeader::HasSequence(const std::string& name) const
{ return d_->sequenceIdLookup_.find(name) != d_->sequenceIdLookup_.cend(); }

inline size_t BamHeader::NumSequences() const
{ return d_->sequences_.size(); }

inline std::string BamHeader::PacBioBamVersion() const
{ return d_->pacbioBamVersion_; }

inline SequenceInfo BamHeader::Sequence(const int32_t id) const
{ return d_->sequences_.at(id); }

inline std::string BamHeader::SequenceLength(const int32_t id) const
{ return Sequence(id).Length(); }

inline std::string BamHeader::SequenceName(const int32_t id) const
{ return Sequence(id).Name(); }

inline std::vector<SequenceInfo> BamHeader::Sequences() const
{ return d_->sequences_; }

inline std::string BamHeader::SortOrder() const
{ return d_->sortOrder_; }

inline BamHeader& BamHeader::SortOrder(std::string order)
{ d_->sortOrder_ = std::move(order); return *this; }

inline std::string BamHeader::Version() const
{ return d_->version_; }

inline BamHeader& BamHeader::Version(std::string version)
{ d_->version_ = std::move(version); return *this; }

} // namespace BAM
} // namespace PacBio

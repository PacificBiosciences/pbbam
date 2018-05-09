// File Description
/// \file SequenceInfo.inl
/// \brief Inline implementations for the SequenceInfo class.
//
// Author: Derek Barnett

#include "pbbam/SequenceInfo.h"

namespace PacBio {
namespace BAM {

inline bool SequenceInfo::operator==(const SequenceInfo& other) const
{
    return assemblyId_ == other.assemblyId_ &&
           checksum_   == other.checksum_   &&
           length_     == other.length_     &&
           name_       == other.name_       &&
           species_    == other.species_    &&
           uri_        == other.uri_        &&
           custom_     == other.custom_;
}

inline bool SequenceInfo::operator!=(const SequenceInfo& other) const
{ return !(*this == other); }

inline std::string SequenceInfo::AssemblyId() const
{ return assemblyId_; }

inline SequenceInfo& SequenceInfo::AssemblyId(std::string id)
{ assemblyId_ = std::move(id); return *this; }

inline std::string SequenceInfo::Checksum() const
{ return checksum_; }

inline SequenceInfo& SequenceInfo::Checksum(std::string checksum)
{ checksum_ = std::move(checksum); return *this; }

inline std::map<std::string, std::string> SequenceInfo::CustomTags() const
{ return custom_; }

inline SequenceInfo& SequenceInfo::CustomTags(std::map<std::string, std::string> custom)
{ custom_ = std::move(custom); return *this; }

inline std::string SequenceInfo::Length() const
{ return length_; }

inline SequenceInfo& SequenceInfo::Length(std::string length)
{ length_ = std::move(length); return *this; }

inline std::string SequenceInfo::Name() const
{ return name_; }

inline SequenceInfo& SequenceInfo::Name(std::string name)
{ name_ = std::move(name); return *this; }

inline std::string SequenceInfo::Species() const
{ return species_; }

inline SequenceInfo& SequenceInfo::Species(std::string species)
{ species_ = std::move(species); return *this; }

inline std::string SequenceInfo::ToSam(const SequenceInfo& seq)
{ return seq.ToSam(); }

inline std::string SequenceInfo::Uri() const
{ return uri_; }

inline SequenceInfo& SequenceInfo::Uri(std::string uri)
{ uri_ = std::move(uri); return *this; }

} // namespace BAM
} // namespace PacBio

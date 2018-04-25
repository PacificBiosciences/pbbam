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

inline SequenceInfo& SequenceInfo::AssemblyId(const std::string& id)
{ assemblyId_ = id; return *this; }

inline std::string SequenceInfo::Checksum() const
{ return checksum_; }

inline SequenceInfo& SequenceInfo::Checksum(const std::string& checksum)
{ checksum_ = checksum; return *this; }

inline std::map<std::string, std::string> SequenceInfo::CustomTags() const
{ return custom_; }

inline SequenceInfo& SequenceInfo::CustomTags(const std::map<std::string, std::string>& custom)
{ custom_ = custom; return *this; }

inline std::string SequenceInfo::Length() const
{ return length_; }

inline SequenceInfo& SequenceInfo::Length(const std::string& length)
{ length_ = length; return *this; }

inline std::string SequenceInfo::Name() const
{ return name_; }

inline SequenceInfo& SequenceInfo::Name(const std::string& name)
{ name_ = name; return *this; }

inline std::string SequenceInfo::Species() const
{ return species_; }

inline SequenceInfo& SequenceInfo::Species(const std::string& species)
{ species_ = species; return *this; }

inline std::string SequenceInfo::ToSam(const SequenceInfo& seq)
{ return seq.ToSam(); }

inline std::string SequenceInfo::Uri() const
{ return uri_; }

inline SequenceInfo& SequenceInfo::Uri(const std::string& uri)
{ uri_ = uri; return *this; }

} // namespace BAM
} // namespace PacBio

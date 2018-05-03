// File Description
/// \file ProgramInfo.inl
/// \brief Inline implementations for the ProgramInfo class.
//
// Author: Derek Barnett

#include "pbbam/ProgramInfo.h"

namespace PacBio {
namespace BAM {

inline std::string ProgramInfo::CommandLine() const
{ return commandLine_; }

inline ProgramInfo& ProgramInfo::CommandLine(std::string cmd)
{ commandLine_ = std::move(cmd); return *this; }

inline std::map<std::string, std::string> ProgramInfo::CustomTags() const
{ return custom_; }

inline ProgramInfo& ProgramInfo::CustomTags(std::map<std::string,
                                            std::string> custom)
{ custom_ = std::move(custom); return *this; }

inline std::string ProgramInfo::Description() const
{ return description_; }

inline ProgramInfo& ProgramInfo::Description(std::string description)
{ description_ = std::move(description); return *this; }

inline std::string ProgramInfo::Id() const
{ return id_; }

inline ProgramInfo& ProgramInfo::Id(std::string id)
{ id_ = std::move(id); return *this; }

inline bool ProgramInfo::IsValid() const
{ return !id_.empty(); }

inline std::string ProgramInfo::Name() const
{ return name_; }

inline ProgramInfo& ProgramInfo::Name(std::string name)
{ name_ = std::move(name); return *this; }

inline std::string ProgramInfo::PreviousProgramId() const
{ return previousProgramId_; }

inline ProgramInfo& ProgramInfo::PreviousProgramId(std::string id)
{ previousProgramId_ = std::move(id); return *this; }

inline std::string ProgramInfo::ToSam(const ProgramInfo& prog)
{ return prog.ToSam(); }

inline std::string ProgramInfo::Version() const
{ return version_; }

inline ProgramInfo& ProgramInfo::Version(std::string version)
{ version_ = std::move(version); return *this; }

} // namespace BAM
} // namespace PacBio

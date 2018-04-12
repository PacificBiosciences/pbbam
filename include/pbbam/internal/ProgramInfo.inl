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

inline ProgramInfo& ProgramInfo::CommandLine(const std::string& cmd)
{ commandLine_ = cmd; return *this; }

inline std::map<std::string, std::string> ProgramInfo::CustomTags() const
{ return custom_; }

inline ProgramInfo& ProgramInfo::CustomTags(const std::map<std::string,
                                            std::string>& custom)
{ custom_ = custom; return *this; }

inline std::string ProgramInfo::Description() const
{ return description_; }

inline ProgramInfo& ProgramInfo::Description(const std::string& description)
{ description_ = description; return *this; }

inline std::string ProgramInfo::Id() const
{ return id_; }

inline ProgramInfo& ProgramInfo::Id(const std::string& id)
{ id_ = id; return *this; }

inline bool ProgramInfo::IsValid() const
{ return !id_.empty(); }

inline std::string ProgramInfo::Name() const
{ return name_; }

inline ProgramInfo& ProgramInfo::Name(const std::string& name)
{ name_ = name; return *this; }

inline std::string ProgramInfo::PreviousProgramId() const
{ return previousProgramId_; }

inline ProgramInfo& ProgramInfo::PreviousProgramId(const std::string& id)
{ previousProgramId_ = id; return *this; }

inline std::string ProgramInfo::ToSam(const ProgramInfo& prog)
{ return prog.ToSam(); }

inline std::string ProgramInfo::Version() const
{ return version_; }

inline ProgramInfo& ProgramInfo::Version(const std::string& version)
{ version_ = version; return *this; }

} // namespace BAM
} // namespace PacBio

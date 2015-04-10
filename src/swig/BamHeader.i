/* BamHeader.i */
%module PacBioBam

%{
#include <pbbam/BamHeader.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::BamHeader::BamHeader(BamHeader&&);      // move ctors not used
%ignore PacBio::BAM::BamHeader::operator=;                   // assignment operators not used

%template(ProgramInfoList)   std::vector<PacBio::BAM::ProgramInfo>;
%template(ReadGroupInfoList) std::vector<PacBio::BAM::ReadGroupInfo>;
%template(SequenceInfoList)  std::vector<PacBio::BAM::SequenceInfo>;

namespace PacBio {
namespace BAM {

class BamHeader
{
public:

#ifdef SWIGR
typedef boost::shared_ptr<PacBio::BAM::BamHeader> SharedPtr;
#else
typedef std::shared_ptr<PacBio::BAM::BamHeader> SharedPtr;
#endif
    
    BamHeader(void);
	BamHeader(const std::string& sam);
    BamHeader(const BamHeader& other);
    ~BamHeader(void);

    std::string PacBioBamVersion(void) const;
    std::string SortOrder(void) const;
    std::string Version(void) const;

    bool HasReadGroup(const std::string& id) const;
    ReadGroupInfo ReadGroup(const std::string& id) const;
    std::vector<std::string> ReadGroupIds(void) const;
    std::vector<ReadGroupInfo> ReadGroups(void) const;

    bool HasSequence(const std::string& name) const;
    int32_t SequenceId(const std::string& name) const;
    std::string SequenceLength(const int32_t id) const;
    std::string SequenceName(const int32_t id) const;
    std::vector<std::string> SequenceNames(void) const;
    SequenceInfo Sequence(const std::string& name) const;
    std::vector<SequenceInfo> Sequences(void) const;

    bool HasProgram(const std::string& id) const;
    ProgramInfo Program(const std::string& id) const;
    std::vector<std::string> ProgramIds(void) const;
    std::vector<ProgramInfo> Programs(void) const;

    std::vector<std::string> Comments(void) const;

    std::string ToSam(void) const;

    void PacBioBamVersion(const std::string& version);
    void SortOrder(const std::string& order);
    void Version(const std::string& version);

    void AddReadGroup(const ReadGroupInfo& readGroup);
    void ClearReadGroups(void);
    void ReadGroups(const std::vector<ReadGroupInfo>& readGroups);

    void AddSequence(const SequenceInfo& sequence);
    void ClearSequences(void);
    void Sequences(const std::vector<SequenceInfo>& sequences);

    void AddProgram(const ProgramInfo& pg);
    void ClearPrograms(void);
    void Programs(const std::vector<ProgramInfo>& programs);

    void AddComment(const std::string& comment);
    void ClearComments(void);
    void Comments(const std::vector<std::string>& comments);
};

} // namespace BAM
} // namespace PacBio

%shared_ptr(PacBio::BAM::BamHeader);

#ifdef SWIGR
%template(BamHeaderSmartPtr) boost::shared_ptr<PacBio::BAM::BamHeader>;
#else
%template(BamHeaderSmartPtr) std::shared_ptr<PacBio::BAM::BamHeader>;
#endif

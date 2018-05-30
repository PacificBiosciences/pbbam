// Author: Derek Barnett

#ifndef PBBAM_VCF_VCFFORMAT_H
#define PBBAM_VCF_VCFFORMAT_H

#include <iosfwd>
#include <string>

#include <pbbam/vcf/VcfHeader.h>
#include <pbbam/vcf/VcfVariant.h>

namespace PacBio {
namespace VCF {

struct VcfFormat
{
public:
    /// \name General format info
    /// \{

    static const char* CurrentVersion();

    /// \}

public:
    /// \name Header format
    /// \{

    static VcfHeader ParsedHeader(const std::string& text);

    static std::string FormattedHeader(const VcfHeader& header);

    static VcfHeader HeaderFromFile(const std::string& fn);

    static VcfHeader HeaderFromStream(std::istream& in);

    /// \}

public:
    /// \name Variant format
    /// \{

    static VcfVariant ParsedVariant(const std::string& line);

    static std::string FormattedVariant(const VcfVariant& var);

    /// \}

    // ---------------------------------------------------------------------- //
    // The following methods are mostly internal helpers, exposed here for    //
    // testing. Client code should probably not need these, but are available //
    // here if needed.                                                        //
    // ---------------------------------------------------------------------- //

public:
    /// \internal
    /// \name Header format helpers
    /// \{

    static ContigDefinition ParsedContigDefinition(std::string line);

    static FilterDefinition ParsedFilterDefinition(std::string line);

    static FormatDefinition ParsedFormatDefinition(std::string line);

    static GeneralDefinition ParsedGeneralDefinition(const std::string& line);

    static InfoDefinition ParsedInfoDefinition(std::string line);

    static std::string FormattedContigDefinition(const ContigDefinition& def);

    static std::string FormattedFilterDefinition(const FilterDefinition& def);

    static std::string FormattedFormatDefinition(const FormatDefinition& def);

    static std::string FormattedGeneralDefinition(const GeneralDefinition& def);

    static std::string FormattedInfoDefinition(const InfoDefinition& def);

    /// \}

public:
    /// \internal
    /// \name Variant format helpers
    /// \{

    static std::string FormattedInfoField(const InfoField& field);

    static std::string FormattedInfoFields(const std::vector<InfoField>& fields);

    static InfoField ParsedInfoField(const std::string& text);

    static std::vector<InfoField> ParsedInfoFields(const std::string& text);

    /// \}
};

}  // namespace VCF
}  // namespace PacBio

#endif  // PBBAM_VCF_VCFFORMAT_H

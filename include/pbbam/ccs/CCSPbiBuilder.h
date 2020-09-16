#ifndef PBBAM_CCS_CCSPBIBUILDER_H
#define PBBAM_CCS_CCSPBIBUILDER_H

#include <pbbam/Config.h>

#include <memory>
#include <string>

#include <pbbam/PbiBuilder.h>
#include <pbbam/ccs/CCSHeader.h>
#include <pbbam/ccs/CCSRecord.h>

namespace PacBio {
namespace CCS {

struct CCSPbiBuilderConfig
{
    using PbiBuilder = BAM::PbiBuilder;

    // zlib compression level for PBI file
    PbiBuilder::CompressionLevel CompressionLevel = PbiBuilder::DefaultCompression;

    // Number of threads to use in PBI file compression. Only active during
    // CCSPbiBuilder::Close().
    size_t NumThreads = 4;
};

class CCSPbiBuilder
{
public:
    CCSPbiBuilder(const std::string& pbiFilename, const std::string& movieName,
                  const CCSPbiBuilderConfig& config = CCSPbiBuilderConfig());
    CCSPbiBuilder(const std::string& pbiFilename, const CCSHeader& header,
                  const CCSPbiBuilderConfig& config = CCSPbiBuilderConfig());
    ~CCSPbiBuilder();

public:
    void AddRecord(const CCSRecord& record);
    void Close();
    const std::string& MovieName() const;

private:
    class CCSPbiBuilderPrivate;
    std::unique_ptr<CCSPbiBuilderPrivate> d_;
};

}  // namespace CCS
}  // namespace PacBio

#endif  //  PBBAM_CCS_CCSPBIBUILDER_H

// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "FileProducer.h"

#include <cstdio>
#include <exception>

namespace PacBio {
namespace BAM {
namespace internal {

FileProducer::FileProducer(std::string targetFilename)
    : FileProducer(std::move(targetFilename), targetFilename + ".tmp")
{
}

FileProducer::FileProducer(std::string targetFilename, std::string tempFilename)
    : targetFilename_{std::move(targetFilename)}, tempFilename_{std::move(tempFilename)}
{
    // override renaming if writing to stdout
    //
    // setting temp filename to '-' keeps consistent interfaces
    // for derived classes to actually operate on temp filename
    if (targetFilename_ == "-") tempFilename_ = "-";
}

FileProducer::~FileProducer()
{
    // skip renaming if there is a 'live' exception
    // or if writing to stdout
    if ((std::current_exception() == nullptr) && (tempFilename_ != "-")) {
        std::rename(tempFilename_.c_str(), targetFilename_.c_str());
    }
}

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

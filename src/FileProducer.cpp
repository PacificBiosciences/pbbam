// Author: Derek Barnett

#include "PbbamInternalConfig.h"

#include "FileProducer.h"

#include <cstdio>

#include <stdexcept>

namespace PacBio {
namespace BAM {

FileProducer::FileProducer(std::string targetFilename)
    : FileProducer(targetFilename, targetFilename + ".tmp")
{
}

FileProducer::FileProducer(std::string targetFilename, std::string tempFilename)
    : targetFilename_{std::move(targetFilename)}, tempFilename_{std::move(tempFilename)}
{
    if (targetFilename_.empty()) {
        throw std::runtime_error{"FileProducer error: cannot write to file with empty name"};
    }

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

}  // namespace BAM
}  // namespace PacBio

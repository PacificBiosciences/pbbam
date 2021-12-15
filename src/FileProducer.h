#ifndef PBBAM_FILEPRODUCER_H
#define PBBAM_FILEPRODUCER_H

#include <pbbam/Config.h>

#include <string>

#include <cstdio>

namespace PacBio {
namespace BAM {

// The FileProducer class provides functionality for working with a temp
// file until successful destruction of a FileProducer-derived class.
//
// Derived classes should be sure to flush/close the temp file, and the
// FileProducer's destructor will ensure that the temp file will be renamed to
// the target filename.
//
// If destruction is triggered by an exception, no renaming will occur.
//
class FileProducer
{
public:
    FileProducer() = delete;

    // Initializes FileProducer with specified target filename. Temp filename is
    // set to target filename plus ".tmp" suffix.
    explicit FileProducer(std::string targetFilename);

    // Initializes FileProducer with specified target filename & explicit temp
    // filename.
    FileProducer(std::string targetFilename, std::string tempFilename);

    // Renames temp file to target filename.
    //
    // Derived classes should ensure that data is flushed and file handle closed
    // before or during their destructor.
    //
    // Remaming will not occur if there is a 'live' exception being thrown.
    //
    virtual ~FileProducer();

    const std::string& TargetFilename() const { return targetFilename_; }
    const std::string& TempFilename() const { return tempFilename_; }

private:
    std::string targetFilename_;
    std::string tempFilename_;
};

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FILEPRODUCER_H

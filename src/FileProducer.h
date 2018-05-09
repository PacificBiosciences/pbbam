// Author: Derek Barnett

#ifndef FILEPRODUCER_H
#define FILEPRODUCER_H

#include <cstdio>
#include <string>

namespace PacBio {
namespace BAM {
namespace internal {

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

protected:
    FileProducer() = delete;

    // Initializes FileProducer with specified target filename. Temp filename is
    // set to target filename plus ".tmp" suffix.
    FileProducer(std::string targetFilename);

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
    ~FileProducer();

protected:
    const std::string& TargetFilename() const { return targetFilename_; }

    const std::string& TempFilename() const { return tempFilename_; }

private:
    std::string targetFilename_;
    std::string tempFilename_;
};

}  // namespace internal
}  // namespace BAM
}  // namespace PacBio

#endif  // FILEPRODUCER_H

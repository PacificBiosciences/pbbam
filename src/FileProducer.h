// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#ifndef FILEPRODUCER_H
#define FILEPRODUCER_H

#include <string>
#include <stdio.h>

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
class FileProducer {

protected:
    FileProducer(void) = delete;

    // Initializes FileProducer with specified target filename. Temp filename is
    // set to target filename plus ".tmp" suffix.
    FileProducer(const std::string& targetFilename);

    // Initializes FileProducer with specified target filename & explicit temp
    // filename.
    FileProducer(const std::string& targetFilename,
                 const std::string& tempFilename);

    // Renames temp file to target filename.
    //
    // Derived classes should ensure that data is flushed and file handle closed
    // before or during their destructor.
    //
    // Remaming will not occur if there is a 'live' exception being thrown.
    //
    ~FileProducer(void);

protected:
    const std::string& TargetFilename(void) const
    { return targetFilename_; }

    const std::string& TempFilename(void) const
    { return tempFilename_; }

private:
    std::string targetFilename_;
    std::string tempFilename_;
};

} // namespace internal
} // namespace BAM
} // namespace PacBio

#endif // FILEPRODUCER_H

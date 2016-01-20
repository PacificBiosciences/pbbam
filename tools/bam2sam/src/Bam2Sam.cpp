// Copyright (c) 2015, Pacific Biosciences of California, Inc.
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
// DISCLAIMED. IN NO EVENT  SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Derek Barnett

#include "Bam2Sam.h"
#include <htslib/sam.h>
#include <stdexcept>
#include <memory>
#include <cassert>
using namespace bam2sam;
using namespace std;

namespace bam2sam {

struct HtslibFileDeleter
{
    void operator()(samFile* file)
    {
        if (file)
            sam_close(file);
        file = nullptr;
    }
};

struct HtslibHeaderDeleter
{
    void operator()(bam_hdr_t* hdr)
    {
        if (hdr)
            bam_hdr_destroy(hdr);
        hdr = nullptr;
    }
};

struct HtslibRecordDeleter
{
    void operator()(bam1_t* b)
    {
        if (b)
            bam_destroy1(b);
        b = nullptr;
    }
};

} // namespace bam2sam

void PbBam2Sam::Run(const Settings &settings)
{
    int htslibResult = 0;

    // open files

    unique_ptr<samFile, HtslibFileDeleter> inFileWrapper(sam_open(settings.inputFilename_.c_str(), "rb"));
    samFile* in = inFileWrapper.get();
    if (!in || !in->fp.bgzf)
        throw std::runtime_error("could not read from stdin");

    unique_ptr<samFile, HtslibFileDeleter> outFileWrapper(sam_open("-", "w"));
    samFile* out = outFileWrapper.get();
    if (!out)
        throw std::runtime_error("could not write to stdout");

    // fetch & write header

    unique_ptr<bam_hdr_t, HtslibHeaderDeleter> headerWrapper(bam_hdr_read(in->fp.bgzf));
    bam_hdr_t* hdr = headerWrapper.get();
    if (!hdr)
        throw std::runtime_error("could not read header");

    if (!settings.noHeader_) {
        htslibResult = sam_hdr_write(out, hdr);
        if (htslibResult != 0)
            throw std::runtime_error("could not write header");
        if (settings.printHeaderOnly_)
            return;
    }

    // fetch & write records

    unique_ptr<bam1_t, HtslibRecordDeleter> recordWrapper(bam_init1());
    bam1_t* b = recordWrapper.get();

    while ((htslibResult = sam_read1(in, hdr, b)) >= 0) {
        htslibResult = sam_write1(out, hdr, b);
        if (htslibResult < 0)
            throw std::runtime_error("error writing record to stdout");
    }
}

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

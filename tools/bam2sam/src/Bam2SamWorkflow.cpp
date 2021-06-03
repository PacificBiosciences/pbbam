#include "Bam2SamWorkflow.h"

#include <memory>
#include <stdexcept>

#include <htslib/sam.h>

#include "Bam2SamSettings.h"

namespace PacBio {
namespace Bam2Sam {
namespace {

struct HtslibFileDeleter
{
    void operator()(samFile* file)
    {
        if (file) {
            sam_close(file);
        }
        file = nullptr;
    }
};

struct HtslibHeaderDeleter
{
    void operator()(bam_hdr_t* hdr)
    {
        if (hdr) {
            bam_hdr_destroy(hdr);
        }
        hdr = nullptr;
    }
};

struct HtslibRecordDeleter
{
    void operator()(bam1_t* b)
    {
        if (b) {
            bam_destroy1(b);
        }
        b = nullptr;
    }
};

}  // namespace

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};

    int htslibResult = 0;

    // open files

    std::unique_ptr<samFile, HtslibFileDeleter> inFileWrapper(
        sam_open(settings.InputFilename.c_str(), "rb"));
    samFile* in = inFileWrapper.get();
    if (!in || !in->fp.bgzf) {
        throw std::runtime_error("could not read from stdin");
    }

    std::unique_ptr<samFile, HtslibFileDeleter> outFileWrapper(sam_open("-", "w"));
    samFile* out = outFileWrapper.get();
    if (!out) {
        throw std::runtime_error("could not write to stdout");
    }

    // fetch & write header

    std::unique_ptr<bam_hdr_t, HtslibHeaderDeleter> headerWrapper(bam_hdr_read(in->fp.bgzf));
    bam_hdr_t* hdr = headerWrapper.get();
    if (!hdr) {
        throw std::runtime_error("could not read header");
    }

    if (!settings.NoHeader) {
        htslibResult = sam_hdr_write(out, hdr);
        if (htslibResult != 0) {
            throw std::runtime_error("could not write header");
        }
        if (settings.HeaderOnly) {
            return EXIT_SUCCESS;
        }
    }

    // fetch & write records

    std::unique_ptr<bam1_t, HtslibRecordDeleter> recordWrapper(bam_init1());
    bam1_t* b = recordWrapper.get();

    while ((htslibResult = sam_read1(in, hdr, b)) >= 0) {
        htslibResult = sam_write1(out, hdr, b);
        if (htslibResult < 0) {
            throw std::runtime_error("error writing record to stdout");
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace Bam2Sam
}  // namespace PacBio

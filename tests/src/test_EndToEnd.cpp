#include <cstdio>
#include <cstdlib>

#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <htslib/sam.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamWriter.h>
#include <pbbam/EntireFileQuery.h>

#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

namespace EndToEndTests {

struct Bam1Deleter
{
    void operator()(bam1_t* b) const noexcept
    {
        if (b != nullptr) {
            bam_destroy1(b);
        }
        b = nullptr;
    }
};

struct SamFileDeleter
{
    void operator()(samFile* file) const noexcept
    {
        if (file != nullptr) {
            sam_close(file);
        }
        file = nullptr;
    }
};

struct BamHdrDeleter
{
    void operator()(bam_hdr_t* hdr) const noexcept
    {
        if (hdr != nullptr) {
            bam_hdr_destroy(hdr);
        }
        hdr = nullptr;
    }
};

const std::string inputBamFn = PbbamTestsConfig::Data_Dir + "/aligned.bam";
const std::string goldStandardSamFn = PbbamTestsConfig::Data_Dir + "/aligned.sam";
const std::string generatedBamFn = PbbamTestsConfig::GeneratedData_Dir + "/generated.bam";
const std::string generatedSamFn = PbbamTestsConfig::GeneratedData_Dir + "/generated.sam";
const std::vector<std::string> generatedFiles = {generatedBamFn, generatedSamFn};

int RunBam2Sam(const std::string& bamFn, const std::string& samFn)
{
    std::ostringstream s;
    s << PbbamTestsConfig::Bam2Sam << " " << bamFn << " > " << samFn;
    return std::system(s.str().c_str());
}

int RunDiff(const std::string& fn1, const std::string& fn2)
{
    std::ostringstream s;
    s << "diff " << fn1 << " " << fn2;
    return std::system(s.str().c_str());
}

void Remove(const std::vector<std::string>& files)
{
    for (const auto& fn : files) {
        std::remove(fn.c_str());
    }
}

void CheckGeneratedOutput()
{
    // convert to sam & diff against gold standard
    const int convertRet = RunBam2Sam(generatedBamFn, generatedSamFn);
    const int diffRet = RunDiff(goldStandardSamFn, generatedSamFn);
    EXPECT_EQ(0, convertRet);
    EXPECT_EQ(0, diffRet);

    // clean up
    Remove(generatedFiles);
}

}  // namespace EndToEndTests

// sanity check for rest of tests below
TEST(BAM_EndToEnd, sanity_check_using_htslib_api_directly)
{
    {  // scoped to force flush & close before conversion/diff

        // open files

        std::unique_ptr<samFile, EndToEndTests::SamFileDeleter> inWrapper{
            sam_open(EndToEndTests::inputBamFn.c_str(), "r")};
        samFile* in = inWrapper.get();
        ASSERT_TRUE(in);

        std::unique_ptr<samFile, EndToEndTests::SamFileDeleter> outWrapper{
            sam_open(EndToEndTests::generatedBamFn.c_str(), "wb")};
        samFile* out = outWrapper.get();
        ASSERT_TRUE(out);

        // fetch & write header

        std::unique_ptr<bam_hdr_t, EndToEndTests::BamHdrDeleter> headerWrapper{sam_hdr_read(in)};
        bam_hdr_t* hdr = headerWrapper.get();
        ASSERT_TRUE(hdr);
        ASSERT_EQ(0, sam_hdr_write(out, hdr));

        // fetch & write records

        std::unique_ptr<bam1_t, EndToEndTests::Bam1Deleter> record{bam_init1()};
        bam1_t* b = record.get();
        ASSERT_TRUE(b);

        while (sam_read1(in, hdr, b) >= 0) {
            const auto ret = sam_write1(out, hdr, b);
            std::ignore = ret;
        }
    }

    EndToEndTests::CheckGeneratedOutput();
}

TEST(BAM_EndToEnd, can_roundtrip_single_thread_count_writing)
{
    {
        const BamFile bamFile{EndToEndTests::inputBamFn};
        BamWriter writer{EndToEndTests::generatedBamFn, bamFile.Header(),
                         BamWriter::DefaultCompression, 1};

        EntireFileQuery entireFile{bamFile};
        for (const BamRecord& record : entireFile) {
            writer.Write(record);
        }
    }

    EndToEndTests::CheckGeneratedOutput();
}

TEST(BAM_EndToEnd, can_roundtrip_default_thread_count_writing)
{
    {
        const BamFile bamFile{EndToEndTests::inputBamFn};
        BamWriter writer{EndToEndTests::generatedBamFn, bamFile.Header()};

        EntireFileQuery entireFile{bamFile};
        for (const BamRecord& record : entireFile) {
            writer.Write(record);
        }
    }

    EndToEndTests::CheckGeneratedOutput();
}

TEST(BAM_EndToEnd, can_roundtrip_system_thread_count_writing)
{
    {
        const BamFile bamFile{EndToEndTests::inputBamFn};
        BamWriter writer{EndToEndTests::generatedBamFn, bamFile.Header(),
                         BamWriter::DefaultCompression, 0};

        EntireFileQuery entireFile{bamFile};
        for (const BamRecord& record : entireFile) {
            writer.Write(record);
        }
    }

    EndToEndTests::CheckGeneratedOutput();
}

TEST(BAM_EndToEnd, can_roundtrip_user_thread_count_writing)
{
    {
        const BamFile bamFile{EndToEndTests::inputBamFn};
        BamWriter writer{EndToEndTests::generatedBamFn, bamFile.Header(),
                         BamWriter::DefaultCompression, 3};

        EntireFileQuery entireFile{bamFile};
        for (const BamRecord& record : entireFile) {
            writer.Write(record);
        }
    }

    EndToEndTests::CheckGeneratedOutput();
}

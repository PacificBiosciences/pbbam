// Author: Derek Barnett

#include <pbbam/FastqWriter.h>

#include <cstddef>
#include <cstdint>

#include <gtest/gtest.h>

#include <pbbam/BamFileMerger.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastqReader.h>
#include <pbbam/FastqSequence.h>

#include "FastxTests.h"
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM_FastqWriter, throws_on_empty_filename)
{
    EXPECT_THROW(FastqWriter writer{""}, std::runtime_error);
}

TEST(BAM_FastqWriter, throws_on_invalid_extension)
{
    EXPECT_THROW(FastqWriter writer{"wrong.ext"}, std::runtime_error);
}

TEST(BAM_FastqWriter, can_write_fastq_sequence)
{
    const std::string outFastq = PbbamTestsConfig::GeneratedData_Dir + "/out.fq";
    const FastqSequence seq{"name", "GATTACA", "!!!!!!!"};

    {
        FastqWriter writer{outFastq};
        writer.Write(seq);
    }

    const auto seqs = FastqReader::ReadAll(outFastq);
    ASSERT_EQ(1, seqs.size());
    EXPECT_EQ(seq.Name(), seqs[0].Name());
    EXPECT_EQ(seq.Bases(), seqs[0].Bases());
    EXPECT_EQ(seq.Qualities(), seqs[0].Qualities());

    remove(outFastq.c_str());
}

TEST(BAM_FastqWriter, can_write_fastq_from_bam)
{
    const std::string fn = PbbamTestsConfig::Data_Dir + "/unmap1.bam";
    const std::string outFastq = PbbamTestsConfig::GeneratedData_Dir + "/out.fq";

    {
        FastqWriter writer{outFastq};
        EntireFileQuery query{fn};
        for (const auto& bam : query)
            writer.Write(bam);
    }
    const auto seqs = FastqReader::ReadAll(outFastq);
    ASSERT_EQ(1, seqs.size());

    const std::string name{"test/1/0_100"};
    const std::string bases{
        "GATCGCACTGAAAATCTGGATATAGAACGTGTGCAAATGATTGTCTCTACCGTTCCGTAAAAATTATTGCTAATTAGCAATGATTTTAAG"
        "CTAATTAGTT"};
    const std::string quals{
        "CCCCCCCCCCCCCCCCCCCACCCCCACCCCCCCCCCCCB;CCCAACCCCCCCCCCCCCD=B9BCABCBCB>BBBC@B<<@BA;BCC?B>"
        "A<<@(?:4==4"};

    EXPECT_EQ(name, seqs[0].Name());
    EXPECT_EQ(bases, seqs[0].Bases());
    EXPECT_EQ(quals, seqs[0].Qualities().Fastq());

    remove(outFastq.c_str());
}

TEST(BAM_FastqWriter, can_write_fastq_from_strings)
{
    const std::string outFastq = PbbamTestsConfig::GeneratedData_Dir + "/out.fq";
    const std::string name = "name";
    const std::string bases = "GATTACA";
    const std::string quals = "!!!!!!!";

    {
        FastqWriter writer{outFastq};
        writer.Write(name, bases, quals);
    }

    const auto seqs = FastqReader::ReadAll(outFastq);
    ASSERT_EQ(1, seqs.size());
    EXPECT_EQ(name, seqs[0].Name());
    EXPECT_EQ(bases, seqs[0].Bases());

    remove(outFastq.c_str());
}

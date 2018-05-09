// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include "PbbamTestData.h"

#include <pbbam/AlignmentPrinter.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/IndexedFastaReader.h>

using namespace PacBio;
using namespace PacBio::BAM;

namespace AlignmentPrinterTests {

const std::string lambdaFasta = PbbamTestsConfig::Data_Dir + "/lambdaNEB.fa";
const std::string singleInsertionBam = PbbamTestsConfig::Data_Dir + "/aligned.bam";

}  // namespace AlignmentPrinterTests

TEST(AlignmentPrinterTest, Print)
{
    IndexedFastaReader r(AlignmentPrinterTests::lambdaFasta);
    AlignmentPrinter pretty(r);

    BamFile bamFile(AlignmentPrinterTests::singleInsertionBam);
    EntireFileQuery bamQuery(bamFile);
    auto it = bamQuery.begin();

    // funky formatting used to format alignments
    auto expected = std::string{
        "Read        : singleInsertion/100/0_49\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 49\n"
        "Concordance : 0.96\n"
        "\n"
        "5210 : GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGG : 5249\n"
        "       \x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| "
        "|\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| "
        "||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "   0 : GGCTGCAG-GTACAGCGGTCAGGAGGCCAATTGATGCCGG :   39\n"
        "\n"
        "5249 : ACTGGCTGAT : 5259\n"
        "       |\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "  39 : ACTGGCTGAT :   49\n"
        "\n"};

    auto record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));

    expected = std::string{
        "Read        : singleInsertion/200/0_49\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 49\n"
        "Concordance : 0.96\n"
        "\n"
        "5210 : GGCTGCAGTGTACAGCGGTCAGGAGGCC-ATTGATGCCGG : 5249\n"
        "       \x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| "
        "|\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||| "
        "||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "   0 : GGCTGCAG-GTACAGCGGTCAGGAGGCCAATTGATGCCGG :   39\n"
        "\n"
        "5249 : ACTGGCTGAT : 5259\n"
        "       |\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||||\n"
        "  39 : ACTGGCTGAT :   49\n"
        "\n"};

    record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));

    expected = std::string{
        "Read        : singleInsertion/100/0_111\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 59\n"
        "Concordance : 0.951\n"
        "\n"
        "9377 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCG : 9417\n"
        "       "
        "|||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||"
        "\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||  |\n"
        "   0 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGA--G :   38\n"
        "\n"
        "9417 : CAGCACGGT-AACAGCGGCAA : 9437\n"
        "       |||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||| "
        "||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||\n"
        "  38 : CAGCACGGTAAACAGCGGCAA :   59\n"
        "\n"};

    record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));

    expected = std::string{
        "Read        : singleInsertion/100/0_111\n"
        "Reference   : lambda_NEB3011\n"
        "\n"
        "Read-length : 59\n"
        "Concordance : 0.951\n"
        "\n"
        "9377 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGATCG : 9417\n"
        "       "
        "|||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||"
        "\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||||||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m|||  |\n"
        "   0 : AAGTCACCAATGTGGGACGTCCGTCGATGGCAGAAGA--G :   38\n"
        "\n"
        "9417 : CAGCACGGT-AACAGCGGCAA : 9437\n"
        "       |||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||| "
        "||||\x1B[1m\x1B[31m|\x1B[0m\x1B[39;49m||||||\n"
        "  38 : CAGCACGGTAAACAGCGGCAA :   59\n"
        "\n"};

    record = *it++;
    EXPECT_EQ(expected, pretty.Print(record, Orientation::GENOMIC));
}

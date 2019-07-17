// Author: Derek Barnett

#ifndef FASTXTESTS_H
#define FASTXTESTS_H

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/FastaSequence.h>
#include <pbbam/FastqSequence.h>

#include "PbbamTestData.h"

using FastaSequence = PacBio::BAM::FastaSequence;
using FastqSequence = PacBio::BAM::FastqSequence;

namespace FastxTests {

// clang-format off

const std::string fastxDataDir = PacBio::BAM::PbbamTestsConfig::Data_Dir + "/fastx/";

const std::string simpleFastaFn        = fastxDataDir + "simple.fa";
const std::string simpleFastaFaiFn     = fastxDataDir + "simple.fa.fai";
const std::string simpleFastaGzipFn    = fastxDataDir + "simple-gzip.fa.gz";
const std::string simpleFastaBgzfFn    = fastxDataDir + "simple-bgzf.fa.gz";
const std::string simpleFastaBgzfGziFn = fastxDataDir + "simple-bgzf.fa.gz.gzi";

const std::string simpleFastqFn        = fastxDataDir + "simple.fq";
const std::string simpleFastqFaiFn     = fastxDataDir + "simple.fq.fai";
const std::string simpleFastqGzipFn    = fastxDataDir + "simple-gzip.fq.gz";
const std::string simpleFastqBgzfFn    = fastxDataDir + "simple-bgzf.fq.gz";
const std::string simpleFastqBgzfGziFn = fastxDataDir + "simple-bgzf.fq.gz.gzi";

const std::string chunkingFastaFn    = fastxDataDir + "chunking.fa";
const std::string chunkingFastaFaiFn = fastxDataDir + "chunking.fa.fai";
const std::string chunkingFastqFn    = fastxDataDir + "chunking.fq";
const std::string chunkingFastqFaiFn = fastxDataDir + "chunking.fq.fai";

const std::vector<FastaSequence> ExpectedFasta {
    FastaSequence{ "seq1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"},
    FastaSequence{ "seq2", "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"},
    FastaSequence{ "seq3", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"},
    FastaSequence{ "seq4", "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"},
    FastaSequence{ "seq5", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"},
    FastaSequence{ "seq6", "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"},
    FastaSequence{ "seq7", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG"},
    FastaSequence{ "seq8", "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"}
};

const std::vector<FastqSequence> ExpectedFastq {
    FastqSequence{
        "seq1",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
      R"(ZABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~)"},
    FastqSequence{
        "seq2",
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
      R"(~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@)"},
    FastqSequence{
        "seq3",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
      R"(!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_)"},
    FastqSequence{
        "seq4",
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
      R"(_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@?>=<;:9876543210/.-,+*)('&%$#"!)"},
    FastqSequence{
        "seq5",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
      R"(;;>@BCEFGHJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~)"},
    FastqSequence{
        "seq6",
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
      R"(~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJHGFECB@>;;)"},
    FastqSequence{
        "seq7",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
      R"(ZABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~)"},
    FastqSequence{
        "seq8",
        "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA",
      R"(~}|{zyxwvutsrqponmlkjihgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@)"},
};



} // namespace FastxTests

#endif  // FASTXTESTS_H
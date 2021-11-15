#include <pbbam/ReadGroupInfo.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>

#include "PbbamTestData.h"

using namespace PacBio::BAM;

namespace ReadGroupHashingTests {

// clang-format off

// -----------------------------------------------
// IDs
//
//  movie name : m54006_200116_134114
//
// CCS:
//                       old : 550216e7
//       old_barcoded (8--8) : c68e726b/8--8
//   old_barcoded (199--199) : 8d2d0124/199--199
//                       new : 550216e7
//       new barcoded (8--8) : 550216e7/8--8
//   new barcoded (199--199) : 550216e7/199--199
//
// SUBREAD:
//                       old : 0388f94c
//       old_barcoded (8--8) : e93f69d9/8--8
//   old_barcoded (199--199) : 9a04acc8/199--199
//                       new : 0388f94c
//       new barcoded (8--8) : 0388f94c/8--8
//   new barcoded (199--199) : 0388f94c/199--199
// -----------------------------------------------

//
// CCS read groups, using legacy RG ID hash
//
static const std::string ccs_no_barcodes_OLD_HASH_RG{
    "@RG\tID:550216e7\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
    "PU:m54006_200116_134114\tPM:SEQUEL"
};
static const std::string ccs_barcode_8_8_OLD_HASH_RG{
    "@RG\tID:c68e726b/8--8\tPL:PACBIO\t"
    "DS:READTYPE=CCS;BINDINGKIT=100-619-300;SEQUENCINGKIT=100-619-400;"
    "BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;BarcodeHash=foo;"
    "BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:8--8"
};
static const std::string ccs_barcode_199_199_OLD_HASH_RG{
    "@RG\tID:8d2d0124/199--199\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:199--199"
};
//
// CCS read groups, using fixed RG ID hash
//
static const std::string ccs_no_barcodes_NEW_HASH_RG{
    "@RG\tID:550216e7\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
    "PU:m54006_200116_134114\tPM:SEQUEL"
};
static const std::string ccs_barcode_8_8_NEW_HASH_RG{
    "@RG\tID:550216e7/8--8\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:8--8"
};
static const std::string ccs_barcode_199_199_NEW_HASH_RG{
    "@RG\tID:550216e7/199--199\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:199--199"
};
//
// subread read groups, using legacy RG ID hash
//
static const std::string subread_no_barcodes_OLD_HASH_RG{
    "@RG\tID:0388f94c\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100\t"
    "PU:m54006_200116_134114\tPM:SEQUEL"
};
static const std::string subread_barcode_8_8_OLD_HASH_RG{
    "@RG\tID:e93f69d9/8--8\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:8--8"
};
static const std::string subread_barcode_199_199_OLD_HASH_RG{
    "@RG\tID:9a04acc8/199--199\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:199--199"
};
//
// subread read groups, using fixed RG ID hash
//
static const std::string subread_no_barcodes_NEW_HASH_RG{
    "@RG\tID:0388f94c\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100"
    "\tPU:m54006_200116_134114\tPM:SEQUEL"
};
static const std::string subread_barcode_8_8_NEW_HASH_RG{
    "@RG\tID:0388f94c/8--8\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:8--8"
};
static const std::string subread_barcode_199_199_NEW_HASH_RG{
    "@RG\tID:0388f94c/199--199\tPL:PACBIO\tDS:READTYPE=SUBREAD;BINDINGKIT=100-619-300;"
    "SEQUENCINGKIT=100-619-400;BASECALLERVERSION=3.0;FRAMERATEHZ=100;BarcodeFile=foo;"
    "BarcodeHash=foo;BarcodeCount=2;BarcodeMode=Symmetric;BarcodeQuality=Score\t"
    "PU:m54006_200116_134114\tPM:SEQUEL\tBC:199--199"
};

static const std::string unrelated_read_group_RG{
    "@RG\tID:ab118ebd\tPL:PACBIO\tDS:READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;"
    "BINDINGKIT=101-490-800;SEQUENCINGKIT=101-490-900;BASECALLERVERSION=5.0.0;"
    "FRAMERATEHZ=100.000000\tPU:m64011_190228_190319\tPM:SEQUELII\tCM:S/P3-C1/5.0-8M"
};

static const std::string Dir = PbbamTestsConfig::Data_Dir + "/read_groups/";

static const std::string ccs_no_barcodes_OLD_HASH_FILE{Dir     + "old_hash_ccs.bam"};                   // 5
static const std::string ccs_barcode_8_8_OLD_HASH_FILE{Dir     + "old_hash_ccs_barcode_8_8.bam"};       // 4
static const std::string ccs_barcode_199_199_OLD_HASH_FILE{Dir + "old_hash_ccs_barcode_199_199.bam"};   // 3
static const std::string ccs_barcodes_mixed_OLD_HASH_FILE{Dir  + "old_hash_ccs_barcodes_mixed.bam"};    // 10

static const std::string ccs_no_barcodes_NEW_HASH_FILE{Dir     + "new_hash_ccs.bam"};                   // 5
static const std::string ccs_barcode_8_8_NEW_HASH_FILE{Dir     + "new_hash_ccs_barcode_8_8.bam"};       // 4
static const std::string ccs_barcode_199_199_NEW_HASH_FILE{Dir + "new_hash_ccs_barcode_199_199.bam"};   // 3
static const std::string ccs_barcodes_mixed_NEW_HASH_FILE{Dir  + "new_hash_ccs_barcodes_mixed.bam"};    // 10

static const std::string subread_no_barcodes_OLD_HASH_FILE{Dir     + "old_hash_subreads.bam"};                  // 3
static const std::string subread_barcode_8_8_OLD_HASH_FILE{Dir     + "old_hash_subreads_barcode_8_8.bam"};      // 2
static const std::string subread_barcode_199_199_OLD_HASH_FILE{Dir + "old_hash_subreads_barcode_199_199.bam"};  // 1
static const std::string subread_barcodes_mixed_OLD_HASH_FILE{Dir  + "old_hash_subreads_barcodes_mixed.bam"};   // 6

static const std::string subread_no_barcodes_NEW_HASH_FILE{Dir     + "new_hash_subreads.bam"};                  // 3
static const std::string subread_barcode_8_8_NEW_HASH_FILE{Dir     + "new_hash_subreads_barcode_8_8.bam"};      // 2
static const std::string subread_barcode_199_199_NEW_HASH_FILE{Dir + "new_hash_subreads_barcode_199_199.bam"};  // 1
static const std::string subread_barcodes_mixed_NEW_HASH_FILE{Dir  + "new_hash_subreads_barcodes_mixed.bam"};   // 6

// clang-format on

void CheckReadGroupFilter(const std::map<std::string, int>& samReadGroupsCounts,
                          const std::string& fn)
{
    SCOPED_TRACE(fn);

    for (const auto& samReadGroupCount : samReadGroupsCounts) {
        SCOPED_TRACE(samReadGroupCount.first);
        const ReadGroupInfo rg = ReadGroupInfo::FromSam(samReadGroupCount.first);
        const PbiReadGroupFilter filter{rg};
        PbiFilterQuery query{filter, fn};
        EXPECT_EQ(samReadGroupCount.second, query.NumReads());
    }
}

}  // namespace ReadGroupHashingTests

// clang-format off

TEST(BAM_ReadGroupHashing, can_filter_old_bam_with_old_barcode_read_hash)
{
    using namespace ReadGroupHashingTests;

    {
        SCOPED_TRACE("file contains barcodes: none");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     5},
            {ccs_barcode_8_8_OLD_HASH_RG,     0},
            {ccs_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_no_barcodes_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     3},
            {subread_barcode_8_8_OLD_HASH_RG,     0},
            {subread_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_no_barcodes_OLD_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     0},
            {ccs_barcode_8_8_OLD_HASH_RG,     4},
            {ccs_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_8_8_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     0},
            {subread_barcode_8_8_OLD_HASH_RG,     2},
            {subread_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_8_8_OLD_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 199--199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     0},
            {ccs_barcode_8_8_OLD_HASH_RG,     0},
            {ccs_barcode_199_199_OLD_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_199_199_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     0},
            {subread_barcode_8_8_OLD_HASH_RG,     0},
            {subread_barcode_199_199_OLD_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_199_199_OLD_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8, 199-199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     0},
            {ccs_barcode_8_8_OLD_HASH_RG,     4},
            {ccs_barcode_199_199_OLD_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcodes_mixed_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     0},
            {subread_barcode_8_8_OLD_HASH_RG,     2},
            {subread_barcode_199_199_OLD_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcodes_mixed_OLD_HASH_FILE);
    }
}

TEST(BAM_ReadGroupHashing, can_filter_old_bam_with_new_barcode_read_hash)
{
    using namespace ReadGroupHashingTests;

    {
        SCOPED_TRACE("file contains barcodes: none");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     5},
            {ccs_barcode_8_8_NEW_HASH_RG,     0},
            {ccs_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_no_barcodes_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     3},
            {subread_barcode_8_8_NEW_HASH_RG,     0},
            {subread_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_no_barcodes_OLD_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     0},
            {ccs_barcode_8_8_NEW_HASH_RG,     4},
            {ccs_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_8_8_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     0},
            {subread_barcode_8_8_NEW_HASH_RG,     2},
            {subread_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_8_8_OLD_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 199--199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     0},
            {ccs_barcode_8_8_NEW_HASH_RG,     0},
            {ccs_barcode_199_199_NEW_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_199_199_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     0},
            {subread_barcode_8_8_NEW_HASH_RG,     0},
            {subread_barcode_199_199_NEW_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_199_199_OLD_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8, 199-199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     0},
            {ccs_barcode_8_8_NEW_HASH_RG,     4},
            {ccs_barcode_199_199_NEW_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcodes_mixed_OLD_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     0},
            {subread_barcode_8_8_NEW_HASH_RG,     2},
            {subread_barcode_199_199_NEW_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcodes_mixed_OLD_HASH_FILE);
    }
}

TEST(BAM_ReadGroupHashing, can_filter_new_bam_with_old_barcode_read_hash)
{
    using namespace ReadGroupHashingTests;

    {
        SCOPED_TRACE("file contains barcodes: none");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     5},
            {ccs_barcode_8_8_OLD_HASH_RG,     0},
            {ccs_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_no_barcodes_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     3},
            {subread_barcode_8_8_OLD_HASH_RG,     0},
            {subread_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_no_barcodes_NEW_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     0},
            {ccs_barcode_8_8_OLD_HASH_RG,     4},
            {ccs_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_8_8_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     0},
            {subread_barcode_8_8_OLD_HASH_RG,     2},
            {subread_barcode_199_199_OLD_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_8_8_NEW_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 199--199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     0},
            {ccs_barcode_8_8_OLD_HASH_RG,     0},
            {ccs_barcode_199_199_OLD_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_199_199_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     0},
            {subread_barcode_8_8_OLD_HASH_RG,     0},
            {subread_barcode_199_199_OLD_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_199_199_NEW_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8, 199-199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_OLD_HASH_RG,     0},
            {ccs_barcode_8_8_OLD_HASH_RG,     4},
            {ccs_barcode_199_199_OLD_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcodes_mixed_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_OLD_HASH_RG,     0},
            {subread_barcode_8_8_OLD_HASH_RG,     2},
            {subread_barcode_199_199_OLD_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcodes_mixed_NEW_HASH_FILE);
    }
}

TEST(BAM_ReadGroupHashing, can_filter_new_bam_with_new_barcode_read_hash)
{
    using namespace ReadGroupHashingTests;

    {
        SCOPED_TRACE("file contains barcodes: none");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     5},
            {ccs_barcode_8_8_NEW_HASH_RG,     0},
            {ccs_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_no_barcodes_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     3},
            {subread_barcode_8_8_NEW_HASH_RG,     0},
            {subread_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_no_barcodes_NEW_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     0},
            {ccs_barcode_8_8_NEW_HASH_RG,     4},
            {ccs_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_8_8_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     0},
            {subread_barcode_8_8_NEW_HASH_RG,     2},
            {subread_barcode_199_199_NEW_HASH_RG, 0},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_8_8_NEW_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 199--199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     0},
            {ccs_barcode_8_8_NEW_HASH_RG,     0},
            {ccs_barcode_199_199_NEW_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcode_199_199_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     0},
            {subread_barcode_8_8_NEW_HASH_RG,     0},
            {subread_barcode_199_199_NEW_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcode_199_199_NEW_HASH_FILE);
    }
    {
        SCOPED_TRACE("file contains barcodes: 8--8, 199-199");

        const std::map<std::string, int> ccsReadGroupCounts{
            {ccs_no_barcodes_NEW_HASH_RG,     0},
            {ccs_barcode_8_8_NEW_HASH_RG,     4},
            {ccs_barcode_199_199_NEW_HASH_RG, 3},
            {unrelated_read_group_RG,         0},
        };
        CheckReadGroupFilter(ccsReadGroupCounts, ccs_barcodes_mixed_NEW_HASH_FILE);

        const std::map<std::string, int> subreadReadGroupCounts{
            {subread_no_barcodes_NEW_HASH_RG,     0},
            {subread_barcode_8_8_NEW_HASH_RG,     2},
            {subread_barcode_199_199_NEW_HASH_RG, 1},
            {unrelated_read_group_RG,             0},
        };
        CheckReadGroupFilter(subreadReadGroupCounts, subread_barcodes_mixed_NEW_HASH_FILE);
    }
}

// clang-format on

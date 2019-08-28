// Author: Derek Barnett

#include <string>

#include <gtest/gtest.h>

#include <pbbam/BgzipWriter.h>
#include <pbbam/FormatUtils.h>
#include "PbbamTestData.h"

using namespace PacBio;
using namespace PacBio::BAM;

TEST(BgzipWriterTest, writes_bgzf_format_to_file)
{
    const std::string fn{PbbamTestsConfig::GeneratedData_Dir + "/bgzf_writer_out.gz"};

    {
        const std::string data{"Simple output data"};
        BgzipWriter writer{fn};
        writer.Write(data.c_str(), data.size());
        writer.Write(data);
    }

    EXPECT_EQ(HtslibCompression::BGZIP, FormatUtils::CompressionType(fn));
}

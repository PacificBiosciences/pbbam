#include "CcsKineticsBystrandifyWorkflow.h"

#include <cassert>
#include <cstdint>

#include <algorithm>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string/predicate.hpp>

#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/SequenceUtils.h>
#include <pbcopper/utility/Ssize.h>

#include <pbbam/BamHeader.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/ProgramInfo.h>

#include "CcsKineticsBystrandifySettings.h"
#include "CcsKineticsBystrandifyVersion.h"

using namespace std::literals::string_literals;

namespace PacBio {
namespace CcsKineticsBystrandify {

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};

    BAM::BamReader inputBamReader{settings.InputFilename};

    // setup our @PG entry to add to header
    BAM::ProgramInfo ccskineticsbystrandifyProgram;
    ccskineticsbystrandifyProgram.Id("ccs-kinetics-bystrandify-"s + CcsKineticsBystrandify::Version)
        .Name("ccs-kinetics-bystrandify")
        .Version(CcsKineticsBystrandify::Version);
    BAM::BamHeader newHeader{inputBamReader.Header().DeepCopy()};
    newHeader.AddProgram(ccskineticsbystrandifyProgram);

    BAM::BamWriter bamWriter{settings.OutputFilename, newHeader};

    for (const auto& read : inputBamReader) {
        const std::string readName = read.FullName();
        PBLOG_VERBOSE << "Processing " << readName;

        if (read.Type() != BAM::RecordType::CCS)
            throw std::runtime_error{"Read '" + readName + "' is of " + BAM::ToString(read.Type()) +
                                     " type, only CCS reads can be converted"};
        if (read.IsMapped())
            throw std::runtime_error{"Read '" + readName +
                                     "' is aligned, only unaligned CCS reads can be converted"};
        if (read.HasPulseWidth())
            throw std::runtime_error{
                "Read '" + readName +
                "' already has 'pw' tag, have you processed this file already?"};
        if (read.HasIPD())
            throw std::runtime_error{
                "Read '" + readName +
                "' already has 'ip' tag, have you processed this file already?"};

        const BAM::BamRecordImpl& readImpl = read.Impl();
        if (!readImpl.HasTag("fn"))
            throw std::runtime_error{"Read '" + readName + "' is missing 'fn' CCS-Kinetics tag"};
        if (!readImpl.HasTag("fp"))
            throw std::runtime_error{"Read '" + readName + "' is missing 'fp' CCS-Kinetics tag"};
        if (!readImpl.HasTag("fi"))
            throw std::runtime_error{"Read '" + readName + "' is missing 'fi' CCS-Kinetics tag"};
        if (!readImpl.HasTag("rn"))
            throw std::runtime_error{"Read '" + readName + "' is missing 'rn' CCS-Kinetics tag"};
        if (!readImpl.HasTag("rp"))
            throw std::runtime_error{"Read '" + readName + "' is missing 'rp' CCS-Kinetics tag"};
        if (!readImpl.HasTag("ri"))
            throw std::runtime_error{"Read '" + readName + "' is missing 'ri' CCS-Kinetics tag"};

        if (boost::ends_with(readName, "/fwd") || boost::ends_with(readName, "/rev"))
            throw std::runtime_error{"Read '" + readName + "' is already by-strandified"};

        // all necessary fields validated, let's create the individual records
        const int32_t holeNumber = read.HoleNumber();
        const auto snr = read.SignalToNoise();
        const auto rq = read.ReadAccuracy();

        std::string seq = read.Sequence();
        Data::QualityValues quals = read.Qualities();
        assert((quals.empty()) || (quals.size() == seq.size()));

        const BAM::ReadGroupInfo rg = read.ReadGroup();
        const BAM::FrameCodec ipdCodec = rg.IpdCodec();
        const BAM::FrameCodec pwCodec = rg.PulseWidthCodec();

        const int32_t fwdPasses = readImpl.TagValue("fn").ToInt32();
        const Data::Frames fwdIPD{ipdCodec == BAM::FrameCodec::V1
                                      ? Data::Frames::Decode(readImpl.TagValue("fi").ToUInt8Array())
                                      : Data::Frames{readImpl.TagValue("fi").ToUInt16Array()}};
        const Data::Frames fwdPW{pwCodec == BAM::FrameCodec::V1
                                     ? Data::Frames::Decode(readImpl.TagValue("fp").ToUInt8Array())
                                     : Data::Frames{readImpl.TagValue("fp").ToUInt16Array()}};
        assert(((fwdPasses == 0) && (fwdIPD.empty())) ||
               ((fwdPasses > 0) && (fwdIPD.size() == seq.size())));
        assert(((fwdPasses == 0) && (fwdPW.empty())) ||
               ((fwdPasses > 0) && (fwdPW.size() == seq.size())));

        const int32_t revPasses = readImpl.TagValue("rn").ToInt32();
        const Data::Frames revIPD{ipdCodec == BAM::FrameCodec::V1
                                      ? Data::Frames::Decode(readImpl.TagValue("ri").ToUInt8Array())
                                      : Data::Frames{readImpl.TagValue("ri").ToUInt16Array()}};
        const Data::Frames revPW{pwCodec == BAM::FrameCodec::V1
                                     ? Data::Frames::Decode(readImpl.TagValue("rp").ToUInt8Array())
                                     : Data::Frames{readImpl.TagValue("rp").ToUInt16Array()}};
        assert(((revPasses == 0) && (revIPD.empty())) ||
               ((revPasses > 0) && (revIPD.size() == seq.size())));
        assert(((revPasses == 0) && (revPW.empty())) ||
               ((revPasses > 0) && (revPW.size() == seq.size())));

        const auto recordWriter = [&newHeader, ipdCodec, pwCodec, holeNumber, &snr, &rq, &rg,
                                   &bamWriter](
            const std::string& newRecordName, const int32_t numPasses, const std::string& sequence,
            const Data::QualityValues& qvs, const Data::Frames& ipd, const Data::Frames& pw) {

            // trim flanking zeroes from IPD/PW vectors (lack of coverage)
            const auto fromStartIt = std::find_if(std::cbegin(ipd), std::cend(ipd),
                                                  [](const uint16_t val) -> bool { return val; });
            const auto fromEndIt = std::find_if(std::crbegin(ipd), std::crend(ipd),
                                                [](const uint16_t val) -> bool { return val; });

            const int32_t beginCutBases = std::distance(std::cbegin(ipd), fromStartIt);
            const int32_t endCutBases = std::distance(std::crbegin(ipd), fromEndIt);

            // can't have an empty sequence
            assert(Utility::Ssize(sequence) - beginCutBases - endCutBases > 0);

            const std::string newSequence(std::cbegin(sequence) + beginCutBases,
                                          std::cend(sequence) - endCutBases);
            const Data::QualityValues newQVs{
                qvs.empty() ? Data::QualityValues{}
                            : Data::QualityValues(std::cbegin(qvs) + beginCutBases,
                                                  std::cend(qvs) - endCutBases)};

            const std::vector<uint16_t> newIpd(std::cbegin(ipd) + beginCutBases,
                                               std::cend(ipd) - endCutBases);
            assert((newIpd.front() != 0) && (newIpd.back() != 0));

            const std::vector<uint16_t> newPW(std::cbegin(pw) + beginCutBases,
                                              std::cend(pw) - endCutBases);

            assert(newQVs.empty() || (newSequence.size() == newQVs.size()));
            assert(newSequence.size() == newIpd.size());
            assert(newSequence.size() == newPW.size());

            if (std::any_of(std::cbegin(newPW), std::cend(newPW),
                            [](const uint16_t val) { return val == 0; })) {
                PBLOG_WARN << "New read '" << newRecordName << "' has '0' PulseWidths, discarding";
                return;
            }

            BAM::BamRecord newRecord{newHeader};
            auto& newRecordImpl = newRecord.Impl();

            // standard CCS defaults
            newRecordImpl.Bin(0)
                .InsertSize(0)
                .MapQuality(255)
                .MatePosition(-1)
                .MateReferenceId(-1)
                .Position(-1)
                .ReferenceId(-1)
                .Flag(0)
                .SetMapped(false);

            BAM::TagCollection tags;
            tags["np"] = numPasses;
            tags["cx"] = static_cast<int32_t>(Data::LocalContextFlags::ADAPTER_BEFORE) |
                         static_cast<int32_t>(Data::LocalContextFlags::ADAPTER_AFTER);

            newRecordImpl.Name(newRecordName)
                .SetSequenceAndQualities(newSequence, newQVs.Fastq())
                .Tags(tags);

            newRecord.IPD(newIpd, ipdCodec)
                .PulseWidth(newPW, pwCodec)
                .HoleNumber(holeNumber)
                .SignalToNoise(snr)
                .ReadAccuracy(rq)
                .ReadGroup(rg);

            bamWriter.Write(newRecord);
        };

        if (fwdPasses >= settings.MinCoverage) {
            recordWriter(readName + "/fwd", fwdPasses, seq, quals, fwdIPD, fwdPW);
        }

        if (revPasses >= settings.MinCoverage) {
            Utility::ReverseComplementCaseSens(seq);
            std::reverse(std::begin(quals), std::end(quals));

            recordWriter(readName + "/rev", revPasses, seq, quals, revIPD, revPW);
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace CcsKineticsBystrandify
}  // namespace PacBio

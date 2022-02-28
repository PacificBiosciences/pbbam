#include "CcsKineticsBystrandifyWorkflow.h"

#include <cassert>
#include <cstdint>

#include <algorithm>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>
#include <pbcopper/utility/SequenceUtils.h>
#include <pbcopper/utility/Ssize.h>

#include <pbbam/BamHeader.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/IndexedBamWriter.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <pbbam/ProgramInfo.h>

#include "CcsKineticsBystrandifySettings.h"
#include "CcsKineticsBystrandifyVersion.h"

namespace PacBio {
namespace CcsKineticsBystrandify {
namespace {

struct StrandifyTask
{
    std::string InputBamFile;
    std::string OutputBamFile;

    BAM::BamHeader NewHeader;
    std::unique_ptr<BAM::BamReader> Reader;
    std::unique_ptr<BAM::IndexedBamWriter> Writer;
};

struct UserIO
{
    UserIO(const Settings& settings)
    {
        const auto CheckInputFile = [](const std::string& fn) {
            if (!Utility::FileExists(fn)) {
                throw std::runtime_error{"Input file does not exist: '" + fn + "' "};
            }
        };

        const auto CheckOutputFile = [](const std::string& fn) {
            if (Utility::FileExists(fn)) {
                PBLOG_WARN << "Overwriting existing output file: " << fn;
            }
        };

        const auto MakeHeaderFrom = [&](const BAM::BamHeader& inputHeader) {
            // add @PG entry to header
            BAM::ProgramInfo ccskineticsbystrandifyProgram;
            ccskineticsbystrandifyProgram
                .Id("ccs-kinetics-bystrandify-" + CcsKineticsBystrandify::Version)
                .Name("ccs-kinetics-bystrandify")
                .Version(CcsKineticsBystrandify::Version)
                .CommandLine(settings.CLI);

            BAM::BamHeader header = inputHeader.DeepCopy();
            header.AddProgram(ccskineticsbystrandifyProgram);
            return header;
        };

        const auto SetupBamIO = [&]() {
            StrandifyTask task;
            task.InputBamFile = settings.InputFilename;
            task.OutputBamFile = settings.OutputFilename;

            CheckInputFile(task.InputBamFile);
            CheckOutputFile(task.OutputBamFile);
            CheckOutputFile(task.OutputBamFile + ".pbi");

            task.Reader = std::make_unique<BAM::BamReader>(task.InputBamFile);
            task.NewHeader = MakeHeaderFrom(task.Reader->Header());
            task.Writer =
                std::make_unique<BAM::IndexedBamWriter>(task.OutputBamFile, task.NewHeader);

            Tasks.push_back(std::move(task));
        };

        const auto ResolveBamPath = [](std::string originalBamPath, std::string datasetPath) {
            assert(!originalBamPath.empty());

            std::string resolvedBamPath;

            // input BAM:  absolute path in input XML
            // output BAM: absolute path in output XML
            if (originalBamPath[0] == '/') {
                resolvedBamPath = originalBamPath;
            }

            // input BAM:  BAM path is relative to input XML
            // output BAM: use same relative-ness to output XML
            else {
                const auto lastSlash = datasetPath.rfind('/');
                if (lastSlash != std::string::npos) {
                    resolvedBamPath = datasetPath.substr(0, lastSlash + 1);
                }
                resolvedBamPath += originalBamPath;
            }

            return resolvedBamPath;
        };

        const auto SetupXmlIO = [&]() {
            IsXml = true;

            InputDatasetFile = settings.InputFilename;
            OutputDatasetFile = settings.OutputFilename;
            CheckInputFile(*InputDatasetFile);
            CheckOutputFile(*OutputDatasetFile);

            const BAM::DataSet dataset{settings.InputFilename};
            const auto filter = BAM::PbiFilter::FromDataSet(dataset);
            assert(dataset.Type() == BAM::DataSet::CONSENSUS_READ);

            const bool isOutputBam = boost::iends_with(*OutputDatasetFile, ".bam");
            const auto& externalResources = dataset.ExternalResources();
            if (isOutputBam && externalResources.Size() != 1) {
                throw std::runtime_error{
                    "Output is BAM. Input XML must only contain 1 input BAM file"};
            }

            for (const BAM::ExternalResource& ext : externalResources) {
                const std::string& bamFilename = ext.ResourceId();
                if (!boost::iends_with(bamFilename, ".bam")) {
                    continue;
                }

                StrandifyTask task;
                task.InputBamFile = ResolveBamPath(bamFilename, *InputDatasetFile);
                task.OutputBamFile = [&]() {
                    if (isOutputBam) {
                        return *OutputDatasetFile;
                    } else {
                        return ResolveBamPath(bamFilename, *OutputDatasetFile);
                    }
                }();
                boost::ireplace_all(task.OutputBamFile, ".bam", ".bystrand.bam");

                CheckInputFile(task.InputBamFile);
                CheckOutputFile(task.OutputBamFile);
                CheckOutputFile(task.OutputBamFile + ".pbi");

                BAM::BamFile bamFile{task.InputBamFile};
                task.NewHeader = MakeHeaderFrom(bamFile.Header());

                task.Reader =
                    filter.IsEmpty()
                        ? std::make_unique<BAM::BamReader>(std::move(bamFile))
                        : std::make_unique<BAM::PbiIndexedBamReader>(filter, std::move(bamFile));

                task.Writer =
                    std::make_unique<BAM::IndexedBamWriter>(task.OutputBamFile, task.NewHeader);

                Tasks.push_back(std::move(task));
            }
        };

        if (boost::iends_with(settings.InputFilename, ".bam")) {
            SetupBamIO();
        } else if (boost::iends_with(settings.InputFilename, ".consensusreadset.xml")) {
            SetupXmlIO();
        } else {
            throw std::runtime_error{
                "Input type is not supported - must be BAM or ConsensusReadSet XML"};
        }
    }

    void WriteXml(int64_t numBases, int64_t numRecords) const
    {
        assert(InputDatasetFile);
        assert(OutputDatasetFile);

        const BAM::DataSet inputDataset{*InputDatasetFile};

        BAM::ConsensusReadSet dataset;
        dataset.Name(inputDataset.Name());
        dataset.Tags(inputDataset.Tags());
        dataset.Filters(inputDataset.Filters());
        dataset.Metadata(inputDataset.Metadata());
        dataset.Metadata().NumRecords(std::to_string(numRecords));
        dataset.Metadata().TotalLength(std::to_string(numBases));

        for (const auto& task : Tasks) {
            BAM::ExternalResource outputBam{"PacBio.ConsensusReadFile.ConsensusReadBamFile",
                                            task.OutputBamFile};
            BAM::FileIndex pbi{"PacBio.Index.PacBioIndex", task.OutputBamFile + ".pbi"};
            outputBam.FileIndices().Add(pbi);
            dataset.ExternalResources().Add(outputBam);
        }

        dataset.Save(*OutputDatasetFile);
    }

    bool IsXml = false;
    std::optional<std::string> InputDatasetFile;
    std::optional<std::string> OutputDatasetFile;

    std::vector<StrandifyTask> Tasks;
};

void Strandify(StrandifyTask& task, const CcsKineticsBystrandify::Settings& settings,
               int64_t& numBases, int64_t& numRecords)
{
    for (const auto& read : *task.Reader) {
        const std::string readName = read.FullName();
        PBLOG_VERBOSE << "Processing " << readName;

        if (read.Type() != BAM::RecordType::CCS) {
            throw std::runtime_error{"Read '" + readName + "' is of " + BAM::ToString(read.Type()) +
                                     " type, only CCS reads can be converted"};
        }
        if (read.IsMapped()) {
            throw std::runtime_error{"Read '" + readName +
                                     "' is aligned, only unaligned CCS reads can be converted"};
        }
        if (read.HasPulseWidth()) {
            throw std::runtime_error{
                "Read '" + readName +
                "' already has 'pw' tag, have you processed this file already?"};
        }
        if (read.HasIPD()) {
            throw std::runtime_error{
                "Read '" + readName +
                "' already has 'ip' tag, have you processed this file already?"};
        }

        const BAM::BamRecordImpl& readImpl = read.Impl();
        if (!readImpl.HasTag("fn")) {
            throw std::runtime_error{"Read '" + readName + "' is missing 'fn' CCS-Kinetics tag"};
        }
        if (!readImpl.HasTag("fp")) {
            throw std::runtime_error{"Read '" + readName + "' is missing 'fp' CCS-Kinetics tag"};
        }
        if (!readImpl.HasTag("fi")) {
            throw std::runtime_error{"Read '" + readName + "' is missing 'fi' CCS-Kinetics tag"};
        }
        if (!readImpl.HasTag("rn")) {
            throw std::runtime_error{"Read '" + readName + "' is missing 'rn' CCS-Kinetics tag"};
        }
        if (!readImpl.HasTag("rp")) {
            throw std::runtime_error{"Read '" + readName + "' is missing 'rp' CCS-Kinetics tag"};
        }
        if (!readImpl.HasTag("ri")) {
            throw std::runtime_error{"Read '" + readName + "' is missing 'ri' CCS-Kinetics tag"};
        }

        if (boost::ends_with(readName, "/fwd") || boost::ends_with(readName, "/rev")) {
            throw std::runtime_error{"Read '" + readName + "' is already by-strandified"};
        }

        // all necessary fields validated, let's create the individual records
        const int32_t holeNumber = read.HoleNumber();
        const auto snr = read.SignalToNoise();
        const auto rq = read.ReadAccuracy();

        std::string seq = read.Sequence();
        Data::QualityValues quals = read.Qualities();
        assert((quals.empty()) || (quals.size() == seq.size()));

        const BAM::ReadGroupInfo rg = read.ReadGroup();
        const Data::FrameCodec ipdCodec = rg.IpdCodec();
        const Data::FrameEncoder ipdEncoder = rg.IpdFrameEncoder();
        const Data::FrameCodec pwCodec = rg.PulseWidthCodec();
        const Data::FrameEncoder pwEncoder = rg.PulseWidthFrameEncoder();

        auto IpdFrames = [&](const std::string& name) -> Data::Frames {
            const auto tag = readImpl.TagValue(name);
            return ipdCodec == Data::FrameCodec::RAW ? Data::Frames{tag.ToUInt16Array()}
                                                     : ipdEncoder.Decode(tag.ToUInt8Array());
        };
        auto PwFrames = [&](const std::string& name) -> Data::Frames {
            const auto tag = readImpl.TagValue(name);
            return pwCodec == Data::FrameCodec::RAW ? Data::Frames{tag.ToUInt16Array()}
                                                    : pwEncoder.Decode(tag.ToUInt8Array());
        };

        const int32_t fwdPasses = readImpl.TagValue("fn").ToInt32();
        const Data::Frames fwdIPD = IpdFrames("fi");
        const Data::Frames fwdPW = PwFrames("fp");
        assert(((fwdPasses == 0) && (fwdIPD.empty())) ||
               ((fwdPasses > 0) && (fwdIPD.size() == seq.size())));
        assert(((fwdPasses == 0) && (fwdPW.empty())) ||
               ((fwdPasses > 0) && (fwdPW.size() == seq.size())));

        const int32_t revPasses = readImpl.TagValue("rn").ToInt32();
        const Data::Frames revIPD = IpdFrames("ri");
        const Data::Frames revPW = PwFrames("rp");
        assert(((revPasses == 0) && (revIPD.empty())) ||
               ((revPasses > 0) && (revIPD.size() == seq.size())));
        assert(((revPasses == 0) && (revPW.empty())) ||
               ((revPasses > 0) && (revPW.size() == seq.size())));

        const auto recordWriter = [&task, ipdCodec, pwCodec, holeNumber, &snr, &rq, &rg, &numBases,
                                   &numRecords](
                                      const std::string& newRecordName, const int32_t numPasses,
                                      const std::string& sequence, const Data::QualityValues& qvs,
                                      const Data::Frames& ipd, const Data::Frames& pw) {
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

            BAM::BamRecord newRecord{task.NewHeader};
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

            task.Writer->Write(newRecord);

            ++numRecords;
            numBases += newRecordImpl.SequenceLength();
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
}

}  // namespace

int Workflow::Runner(const CLI_v2::Results& args)
{
    const Settings settings{args};
    UserIO uio{settings};

    int64_t numBases = 0;
    int64_t numRecords = 0;
    for (auto& task : uio.Tasks) {
        Strandify(task, settings, numBases, numRecords);
    }

    if (uio.IsXml) {
        uio.WriteXml(numBases, numRecords);
    }

    return EXIT_SUCCESS;
}

}  // namespace CcsKineticsBystrandify
}  // namespace PacBio

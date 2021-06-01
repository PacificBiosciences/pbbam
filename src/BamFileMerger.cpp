#include "PbbamInternalConfig.h"

#include <pbbam/BamFileMerger.h>

#include <memory>
#include <stdexcept>

#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/CompositeBamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/IRecordWriter.h>
#include <pbbam/IndexedBamWriter.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiIndexedBamReader.h>

namespace PacBio {
namespace BAM {
namespace {

std::unique_ptr<IRecordWriter> MakeBamWriter(BamHeader header, const std::string& outputFilename,
                                             const bool createPbi, const ProgramInfo& pgInfo)
{
    if (outputFilename.empty()) {
        throw std::runtime_error{"[pbbam] BAM file merging ERROR: no output filename provided"};
    }

    if (pgInfo.IsValid()) {
        header.AddProgram(pgInfo);
    }

    // make BAM writer
    if (createPbi) {
        return std::make_unique<IndexedBamWriter>(outputFilename, header);
    } else {
        return std::make_unique<BamWriter>(outputFilename, header);
    }
}

template <typename Reader>
void MergeImpl(IRecordWriter& writer, Reader& reader)
{
    BamRecord record;
    while (reader.GetNext(record)) {
        writer.Write(record);
    }
}

void MergeToWriter(const DataSet& dataset, const BamHeader& header, IRecordWriter& writer)
{
    const bool isCoordinateSorted = (header.SortOrder() == "coordinate");
    const auto filter = PbiFilter::FromDataSet(dataset);

    if (isCoordinateSorted) {
        if (filter.IsEmpty()) {
            SortedCompositeBamReader<Compare::AlignmentPosition> reader{dataset};
            MergeImpl(writer, reader);
        } else {
            PbiFilterCompositeBamReader<Compare::AlignmentPosition> reader{filter, dataset};
            MergeImpl(writer, reader);
        }
    } else {
        if (filter.IsEmpty()) {
            SortedCompositeBamReader<Compare::QName> reader{dataset};
            MergeImpl(writer, reader);
        } else {
            PbiFilterCompositeBamReader<Compare::QName> reader{filter, dataset};
            MergeImpl(writer, reader);
        }
    }
}

void MergeToWriter(const DataSet& dataset, IRecordWriter& writer)
{
    MergeToWriter(dataset, BamHeader{dataset}, writer);
}

void MergeToFile(const DataSet& dataset, const std::string& outputFilename, bool createPbi,
                 const ProgramInfo& pgInfo)
{
    const BamHeader header{dataset};
    auto writer = MakeBamWriter(header, outputFilename, createPbi, pgInfo);
    MergeToWriter(dataset, header, *writer);
}

}  // namespace

void BamFileMerger::Merge(const DataSet& dataset, const std::string& outputFilename, bool createPbi,
                          const ProgramInfo& pgInfo)
{
    MergeToFile(dataset, outputFilename, createPbi, pgInfo);
}

void BamFileMerger::Merge(const std::vector<std::string>& bamFilenames,
                          const std::string& outputFilename, bool createPbi,
                          const ProgramInfo& pgInfo)
{
    Merge(DataSet{bamFilenames}, outputFilename, createPbi, pgInfo);
}

void BamFileMerger::Merge(const DataSet& dataset, IRecordWriter& writer)
{
    MergeToWriter(dataset, writer);
}

void BamFileMerger::Merge(const std::vector<std::string>& bamFilenames, IRecordWriter& writer)
{
    Merge(DataSet{bamFilenames}, writer);
}

}  // namespace BAM
}  // namespace PacBio

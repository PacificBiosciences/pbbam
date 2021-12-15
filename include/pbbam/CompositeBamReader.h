#ifndef PBBAM_COMPOSITEBAMREADER_H
#define PBBAM_COMPOSITEBAMREADER_H

#include <pbbam/Config.h>

#include <pbbam/BaiIndexCache.h>
#include <pbbam/BaiIndexedBamReader.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/BamReader.h>
#include <pbbam/BamRecord.h>
#include <pbbam/DataSet.h>
#include <pbbam/GenomicInterval.h>
#include <pbbam/PbiIndexedBamReader.h>

#include <deque>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

namespace internal {

/// \internal
/// \brief The CompositeMergeItem class provides a helper struct for composite
///        readers, containing a single-file reader and its "next" record.
///
struct CompositeMergeItem
{
public:
    std::unique_ptr<BamReader> reader;
    BamRecord record;

public:
    CompositeMergeItem(std::unique_ptr<BamReader> rdr);
    CompositeMergeItem(std::unique_ptr<BamReader> rdr, BamRecord rec);
};

/// \internal
/// \brief The CompositeMergeItemSorter class provides a helper function object
///        for ordering composite reader results.
///
/// Essentially just exracts a BamRecord from its parent CompositeMergeItem for
/// further checks.
///
template <typename CompareType>
struct CompositeMergeItemSorter
{
    bool operator()(const CompositeMergeItem& lhs, const CompositeMergeItem& rhs) const;
};

}  // namespace internal

struct PositionSorter  //: public CompositeMergeItemSorter<Compare::AlignmentPosition>
{
    bool operator()(const internal::CompositeMergeItem& lhs,
                    const internal::CompositeMergeItem& rhs) const
    {
        return cmp_(lhs.record, rhs.record);
    }

    Compare::AlignmentPosition cmp_;
};

struct QNameSorter
{
    bool operator()(const internal::CompositeMergeItem& lhs,
                    const internal::CompositeMergeItem& rhs) const
    {
        return cmp_(lhs.record, rhs.record);
    }

    Compare::QName cmp_;
};

template <typename OrderByType>
class PBBAM_EXPORT SortedCompositeBamReader : public internal::IQuery
{
public:
    using value_type = internal::CompositeMergeItem;
    using merge_sorter_type = internal::CompositeMergeItemSorter<OrderByType>;
    using container_type = std::multiset<value_type, merge_sorter_type>;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

public:
    SortedCompositeBamReader(const DataSet& dataset);
    SortedCompositeBamReader(std::vector<BamFile> bamFiles);

    SortedCompositeBamReader(SortedCompositeBamReader&&) noexcept;
    SortedCompositeBamReader& operator=(SortedCompositeBamReader&&) noexcept;
    ~SortedCompositeBamReader() override;

    bool GetNext(BamRecord& record) override;

protected:
    std::vector<BamFile> bamFiles_;
    container_type mergeItems_;  //mergeItems_;
};

/// \brief The GenomicIntervalCompositeBamReader class provides read access to
///        multipe %BAM files, limiting results to a genomic region.
///
/// Requires a ".bai" file for each input %BAM file.
///
/// Results will be returned in order of genomic coordinate (first by reference
/// ID, then by position).
///
class PBBAM_EXPORT GenomicIntervalCompositeBamReader
    : public SortedCompositeBamReader<Compare::AlignmentPosition>
{
public:
    /// \name Contstructors & Related Methods
    /// \{

    /// \brief Constructs composite %BAM reader, that can be queried on genomic interval.
    ///
    /// \note Using this constructor means that an interval must be provided, via
    ///       reader.Interval(i), before iterating.
    ///
    /// \param[in] bamFiles   input BamFile objects
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    GenomicIntervalCompositeBamReader(const std::vector<BamFile>& bamFiles);
    GenomicIntervalCompositeBamReader(const std::vector<BamFile>& bamFiles,
                                      const BaiIndexCache& cache);

    /// \brief Constructs composite %BAM reader, that can be queried on genomic interval.
    ///
    /// \note Using this constructor means that an interval must be provided, via
    ///       reader.Interval(i), before iterating.
    ///
    /// \param[in] dataset      input DataSet
    ///
    /// \throws std::runtime_error if either file (*.bam or *.bai) fails to open
    ///         for reading, or if the interval is invalid
    ///
    GenomicIntervalCompositeBamReader(const DataSet& dataset);
    GenomicIntervalCompositeBamReader(const DataSet& dataset, const BaiIndexCache& cache);

    /// \brief Constructs composite %BAM reader, limiting record results to
    ///        only those overalpping a GenomicInterval.
    ///
    /// \param[in] interval   genomic interval of interest
    /// \param[in] bamFiles   input BamFile objects
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         BAI files.
    ///
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                      const std::vector<BamFile>& bamFiles);
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval,
                                      const std::vector<BamFile>& bamFiles,
                                      const BaiIndexCache& cache);

    /// \brief Constructs composite %BAM reader, limiting record results to
    ///        only those overalpping a GenomicInterval.
    ///
    /// \param[in] interval genomic interval of interest
    /// \param[in] dataset  input data source(s)
    ///
    /// \throws std::runtime_error on failure to open/read underlying %BAM or
    ///         BAI files.
    ///
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval, const DataSet& dataset);
    GenomicIntervalCompositeBamReader(const GenomicInterval& interval, const DataSet& dataset,
                                      const BaiIndexCache& cache);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Sets a new genomic interval of interest.
    ///
    /// \returns reference to this reader
    ///
    GenomicIntervalCompositeBamReader& Interval(const GenomicInterval& interval);

    /// \returns the current specified interval
    ///
    const GenomicInterval& Interval() const;

    /// \}

private:
    BaiIndexCache indexCache_;
    GenomicInterval interval_;
};

/// \brief Provides read access to multipe %BAM files, limiting results to those
///        passing a PbiFilter.
///
/// Requires a ".pbi" file for each input %BAM file.
///
/// \note The template parameter OrderByType is not fully implemented at this
///       time. Use of comparison functor (e.g. Compare::Zmw) for this will
///       currently result in the proper "next" value <b> at each iteration
///       step, independently, but not over the full data set. </b> If all
///       files' "order-by" data values are accessible in increasing order
///       within each file, then the expected ordering will be observed,
///       However, if these data are not sorted within a file, the final results
///       will appear unordered. \n
///       \n
///           Example:\n
///           file 1: { 1, 5, 2, 6 } \n
///           file 2: { 3, 8, 4, 7 } \n
///           results: { 1, 3, 5, 2, 6, 8, 4, 7 } \n
///       \n
///       This a known issue and will be addressed in a future update. But in
///       the meantime, use of Compare::None as the OrderByType is recommended,
///       to explicitly indicate that no particular ordering is expected.
///
template <typename OrderByType = Compare::None>
class PBBAM_EXPORT PbiFilterCompositeBamReader : public SortedCompositeBamReader<OrderByType>
{
public:
    using value_type = internal::CompositeMergeItem;
    using merge_sorter_type = internal::CompositeMergeItemSorter<OrderByType>;
    using container_type = std::multiset<value_type, merge_sorter_type>;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;

public:
    /// \name Contstructors & Related Methods
    /// \{

    PbiFilterCompositeBamReader(const PbiFilter& filter, const std::vector<BamFile>& bamFiles);
    PbiFilterCompositeBamReader(const PbiFilter& filter, const std::vector<BamFile>& bamFiles,
                                const PbiIndexCache& cache);

    PbiFilterCompositeBamReader(const PbiFilter& filter, const DataSet& dataset);
    PbiFilterCompositeBamReader(const PbiFilter& filter, const DataSet& dataset,
                                const PbiIndexCache& cache);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Sets a new PBI filter
    ///
    /// \returns reference to this reader
    ///
    PbiFilterCompositeBamReader& Filter(const PbiFilter& filter);

    /// \returns number of reads that pass the current filter
    ///
    uint32_t NumReads() const;

    /// \}

private:
    PbiIndexCache indexCache_;
    uint32_t numReads_;
};

/// \brief The SequentialCompositeBamReader class provides read access to
///        multiple %BAM files, reading through the entire contents of each
///        file.
///
/// Input files will be accessed in the order provided to the constructor. Each
/// file's contents will be exhausted before moving on to the next one (as
/// opposed to a "round-robin" scheme).
///
class PBBAM_EXPORT SequentialCompositeBamReader : public internal::IQuery
{
public:
    /// \name Contstructors & Related Methods
    /// \{

    SequentialCompositeBamReader(std::vector<BamFile> bamFiles);
    SequentialCompositeBamReader(const DataSet& dataset);

    /// \}

public:
    /// \name Data Access
    /// \{

    /// Fetches next BAM record from the .
    ///
    /// \returns true on success, false if no more data available.
    ///
    bool GetNext(BamRecord& record);

    /// \}

private:
    std::deque<std::unique_ptr<BamReader>> readers_;
};

}  // namespace BAM
}  // namespace PacBio

#include <pbbam/internal/CompositeBamReader.inl>

#endif  // PBBAM_COMPOSITEBAMREADER_H

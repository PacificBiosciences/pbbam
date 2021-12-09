#ifndef PBBAM_FAIINDEX_H
#define PBBAM_FAIINDEX_H

#include <pbbam/Config.h>

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace BAM {

struct FaiEntry
{
    /// Total length of this reference sequence, in bases
    uint64_t Length = 0;

    /// Offset in the FASTA/FASTQ file of this sequence's first base
    uint64_t SeqOffset = 0;

    /// The number of bases on each line
    uint16_t NumBases = 0;

    // The number of bytes in each line, including the newline (allows for Windows newlines)
    uint16_t NumBytes = 0;

    // Offset of sequence's first quality within the FASTQ file (-1 if FASTA only)
    int64_t QualOffset = -1;
};

class FaiIndex
{
public:
    ///
    /// Create *.fai for a FASTA or FASTQ file.
    ///
    static void Create(const std::string& fn);

public:
    ///
    /// \brief Load FAI data from \p fn (*.fai)
    ///
    explicit FaiIndex(const std::string& fn);

    FaiIndex();
    FaiIndex(FaiIndex&&) noexcept;
    FaiIndex& operator=(FaiIndex&&) noexcept;
    ~FaiIndex();

    ///
    /// \brief Add new FAI line to index
    ///
    void Add(std::string name, FaiEntry entry);

    ///
    /// \returns FAI entry for sequence name
    ///
    const FaiEntry& Entry(const std::string& name) const;

    ///
    /// \returns FAI entry at \p row
    ///
    const FaiEntry& Entry(uint32_t row) const;

    ///
    /// \returns true if sequence name found in index
    ///
    bool HasEntry(const std::string& name) const;

    ///
    /// \returns sequence names in FAI file
    ///
    const std::vector<std::string>& Names() const;

    ///
    /// \brief Save FAI data to file
    ///
    void Save(const std::string& fn) const;

    ///
    /// \brief Save FAI data to output stream
    ///
    void Save(std::ostream& out) const;

private:
    class FaiIndexPrivate;
    std::unique_ptr<FaiIndexPrivate> d_;
};

bool operator==(const FaiEntry& lhs, const FaiEntry& rhs) noexcept;

// NOTE: FaiEntry output *does not* include the name column, FaiIndex::Save()
//       handles this mapping
std::ostream& operator<<(std::ostream& out, const FaiEntry& entry);

}  // namespace BAM
}  // namespace PacBio

#endif  // PBBAM_FAIINDEX_H

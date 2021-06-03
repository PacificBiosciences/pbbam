#include "PbbamInternalConfig.h"

#include <pbbam/FastaReader.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <pbbam/FormatUtils.h>

#include "KSeqReader.h"

namespace PacBio {
namespace BAM {

class FastaReader::FastaReaderPrivate
{
public:
    explicit FastaReaderPrivate(const std::string& fn)
    {
        // validate extension
        if (!FormatUtils::IsFastaFilename(fn)) {
            std::ostringstream s;
            s << "[pbbam] FASTA reader ERROR: not a recognized FASTA extension:\n"
              << "  file: " << fn;
            throw std::runtime_error{s.str()};
        }
        reader_ = std::make_unique<KSeqReader>(fn);
    }

    bool GetNext(FastaSequence& record)
    {
        const auto readOk = reader_->ReadNext();
        if (!readOk) {
            return false;  // not error, could be EOF
        }

        record = FastaSequence{reader_->Name(), reader_->Bases()};
        return true;
    }

    std::unique_ptr<KSeqReader> reader_;
};

FastaReader::FastaReader(const std::string& fn)
    : internal::QueryBase<FastaSequence>{}, d_{std::make_unique<FastaReaderPrivate>(fn)}
{
}

FastaReader::FastaReader(FastaReader&&) noexcept = default;

FastaReader& FastaReader::operator=(FastaReader&&) noexcept = default;

FastaReader::~FastaReader() = default;

bool FastaReader::GetNext(FastaSequence& record) { return d_->GetNext(record); }

std::vector<FastaSequence> FastaReader::ReadAll(const std::string& fn)
{
    std::vector<FastaSequence> result;
    result.reserve(256);
    FastaReader reader{fn};
    for (const auto& seq : reader) {
        result.emplace_back(seq);
    }
    return result;
}

}  // namespace BAM
}  // namespace PacBio

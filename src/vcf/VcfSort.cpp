// Author: Derek Barnett

#include <pbbam/vcf/VcfSort.h>

#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pbbam/vcf/VcfQuery.h>
#include <pbbam/vcf/VcfWriter.h>

namespace PacBio {
namespace VCF {

void SortFile(const VcfFile& file, const std::string& outputFilename)
{
    const auto& header = file.Header();

    // configure contig sort order
    std::unordered_map<std::string, size_t> contigLookup;
    const auto& contigDefs = header.ContigDefinitions();
    for (size_t i = 0; i < contigDefs.size(); ++i) {
        const auto& contigId = contigDefs.at(i).Id();
        contigLookup.insert(std::make_pair(contigId, i));
    }

    // read & sort variants
    std::vector<VcfVariant> variants;
    VcfQuery query{file};
    for (const auto& v : query)
        variants.push_back(v);

    std::sort(variants.begin(), variants.end(),
              [&contigLookup](const VcfVariant& lhs, const VcfVariant& rhs) {
                  const auto lhsIdx = contigLookup.at(lhs.Chrom());
                  const auto rhsIdx = contigLookup.at(rhs.Chrom());
                  const auto lhsPos = lhs.Position();
                  const auto rhsPos = rhs.Position();
                  return std::tie(lhsIdx, lhsPos) < std::tie(rhsIdx, rhsPos);
              });

    // write results to file
    VcfWriter writer{outputFilename, header};
    for (const auto& var : variants)
        writer.Write(var);
}

}  // namespace VCF
}  // namespace PacBio

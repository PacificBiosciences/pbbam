#include "PbbamInternalConfig.h"

#include <pbbam/AlignmentPrinter.h>

#include <pbbam/BamRecord.h>
#include <pbbam/IndexedFastaReader.h>

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include <cmath>
#include <cstddef>

namespace PacBio {
namespace BAM {

AlignmentPrinter::AlignmentPrinter(const IndexedFastaReader& ifr)
    : ifr_{std::make_unique<IndexedFastaReader>(ifr)}
{}

std::string AlignmentPrinter::Print(const BamRecord& record, const Data::Orientation orientation)
{
    const std::string seq{record.Sequence(orientation, true, true)};
    const std::string ref{ifr_->ReferenceSubsequence(record, orientation, true, true)};

    if (seq.size() != ref.size()) {
        std::ostringstream s;
        s << "[pbbam] alignment printer ERROR: sequence and reference lengths are not equal:\n"
          << "  seq: " << seq.size() << '\n'
          << "  ref: " << ref.size();
        throw std::runtime_error{s.str()};
    }

    int seqLength = 0;
    float matches = 0;
    std::string pretty;
    Data::Position refCoord = record.ReferenceStart();
    Data::Position seqCoord = BAM::IsCcsOrTranscript(record.Type()) ? 0 : record.QueryStart();

    for (size_t i = 0; i < seq.size();) {
        auto refCoordStr = std::to_string(refCoord);
        auto seqCoordStr = std::to_string(seqCoord);

        size_t maxCoordLength = std::max(refCoordStr.size(), seqCoordStr.size());
        while (refCoordStr.size() < maxCoordLength) {
            refCoordStr = " " + refCoordStr;
        }
        while (seqCoordStr.size() < maxCoordLength) {
            seqCoordStr = " " + seqCoordStr;
        }

        std::string seqWrap{seqCoordStr + " : "};
        std::string refWrap{refCoordStr + " : "};
        std::string prettyWrap(maxCoordLength + 3, ' ');
        prettyWrap.reserve(seq.size());

        // clang-format off
        for (int j = 0; i < seq.size() && j < 40; ++i, ++j) {
            refWrap += ref[i];

            if (seq[i] == ref[i]) {
                ++matches;
                if (refCoord == 0 || refCoord % 10) {
                    prettyWrap += '|';
                } else {
                    prettyWrap += "\033" "[1m" "\x1b" "[31m";
                    prettyWrap += '|';
                    prettyWrap += "\033" "[0m" "\x1b" "[39;49m";
                }
                seqWrap += seq[i];
            } else if (seq[i] == '-' || ref[i] == '-') {
                prettyWrap += ' ';
                seqWrap += seq[i];
            } else {
                prettyWrap += '.';
                seqWrap += "\033" "[1m" "\x1b" "[31m";
                seqWrap += seq[i];
                seqWrap += "\033" "[0m" "\x1b" "[39;49m";
            }
            if (seq[i] != '-') {
                ++seqLength;
                ++seqCoord;
            }
            if (ref[i] != '-') {
                ++refCoord;
            }
        }
        // clang-format on

        refCoordStr = std::to_string(refCoord);
        seqCoordStr = std::to_string(seqCoord);

        maxCoordLength = std::max(refCoordStr.size(), seqCoordStr.size());
        while (refCoordStr.size() < maxCoordLength) {
            refCoordStr = " " + refCoordStr;
        }
        while (seqCoordStr.size() < maxCoordLength) {
            seqCoordStr = " " + seqCoordStr;
        }

        seqWrap += " : " + seqCoordStr;
        refWrap += " : " + refCoordStr;

        pretty += refWrap + '\n' + prettyWrap + '\n' + seqWrap + "\n\n";
    }
    const float similarity = matches / seq.size();

    std::ostringstream output;
    output << "Read        : " << record.FullName() << '\n'
           << "Reference   : " << record.ReferenceName() << "\n\n"
           << "Read-length : " << seqLength << '\n'
           << "Concordance : " << std::setprecision(3) << (similarity) << "\n\n"
           << pretty;
    return output.str();
}

}  // namespace BAM
}  // namespace PacBio

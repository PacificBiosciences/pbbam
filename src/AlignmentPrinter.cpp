// File Description
/// \file AlignmentPrinter.cpp
/// \brief Implements the AlignmentPrinter class.
//
// Author: Armin Töpfer

#include "PbbamInternalConfig.h"

#include "pbbam/AlignmentPrinter.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "pbbam/MakeUnique.h"

namespace PacBio {
namespace BAM {

AlignmentPrinter::AlignmentPrinter(const IndexedFastaReader& ifr)
    : ifr_{std::make_unique<IndexedFastaReader>(ifr)}
{
}

std::string AlignmentPrinter::Print(const BamRecord& record, const Orientation orientation)
{
    const std::string seq{record.Sequence(orientation, true, true)};
    const std::string ref{ifr_->ReferenceSubsequence(record, orientation, true, true)};

    if (seq.size() != ref.size())
        throw std::runtime_error{"Sequence and reference parts are of different size"};

    int seqLength = 0;
    float matches = 0;
    std::string pretty;
    Position refCoord = record.ReferenceStart();
    Position seqCoord = record.QueryStart();

    for (size_t i = 0; i < seq.size();) {
        auto refCoordStr = std::to_string(refCoord);
        auto seqCoordStr = std::to_string(seqCoord);

        size_t maxCoordLength = std::max(refCoordStr.size(), seqCoordStr.size());
        while (refCoordStr.size() < maxCoordLength)
            refCoordStr = " " + refCoordStr;
        while (seqCoordStr.size() < maxCoordLength)
            seqCoordStr = " " + seqCoordStr;

        std::string seqWrap{seqCoordStr + " : "};
        std::string refWrap{refCoordStr + " : "};
        std::string prettyWrap(maxCoordLength + 3, ' ');
        prettyWrap.reserve(seq.size());

        for (int j = 0; i < seq.size() && j < 40; ++i, ++j) {
            refWrap += ref[i];

            if (seq[i] == ref[i]) {
                ++matches;
                if (refCoord == 0 || refCoord % 10)
                    prettyWrap += '|';
                else {
                    prettyWrap += "\033[1m\x1b[31m";
                    prettyWrap += '|';
                    prettyWrap += "\033[0m\x1b[39;49m";
                }
                seqWrap += seq[i];
            } else if (seq[i] == '-' || ref[i] == '-') {
                prettyWrap += ' ';
                seqWrap += seq[i];
            } else {
                prettyWrap += '.';
                seqWrap += "\033[1m\x1b[31m";
                seqWrap += seq[i];
                seqWrap += "\033[0m\x1b[39;49m";
            }
            if (seq[i] != '-') {
                ++seqLength;
                ++seqCoord;
            }
            if (ref[i] != '-') {
                ++refCoord;
            }
        }

        refCoordStr = std::to_string(refCoord);
        seqCoordStr = std::to_string(seqCoord);

        maxCoordLength = std::max(refCoordStr.size(), seqCoordStr.size());
        while (refCoordStr.size() < maxCoordLength)
            refCoordStr = " " + refCoordStr;
        while (seqCoordStr.size() < maxCoordLength)
            seqCoordStr = " " + seqCoordStr;

        seqWrap += " : " + seqCoordStr;
        refWrap += " : " + refCoordStr;

        pretty += refWrap + '\n' + prettyWrap + '\n' + seqWrap + "\n\n";
    }
    const float similarity = matches / seq.size();

    std::stringstream output;

    output << "Read        : " << record.FullName() << std::endl;
    output << "Reference   : " << record.ReferenceName() << std::endl;
    output << std::endl;
    output << "Read-length : " << seqLength << std::endl;
    output << "Concordance : " << std::setprecision(3) << (similarity);
    output << std::endl;
    output << std::endl;
    output << pretty;

    return output.str();
}

}  // namespace BAM
}  // namespace PacBio

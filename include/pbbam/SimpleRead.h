// File Description
/// \file SImpleRead.h
/// \brief Defines the SimpleRead class.
//
// Author: Derek Barnett

#ifndef SIMPLEREAD_H
#define SIMPLEREAD_H

#include "pbbam/Config.h"

#include <string>

#include <boost/optional.hpp>

#include "pbbam/BamRecord.h"
#include "pbbam/Cigar.h"
#include "pbbam/Frames.h"
#include "pbbam/Position.h"
#include "pbbam/QualityValues.h"
#include "pbbam/SNR.h"
#include "pbbam/Strand.h"

namespace PacBio {
namespace BAM {

class SimpleRead
{
public:
    explicit SimpleRead(const BamRecord& bam);
    SimpleRead(std::string name, std::string seq, QualityValues qualities, SNR snr);
    SimpleRead(std::string name, std::string seq, QualityValues qualities, SNR snr, Position qStart,
               Position qEnd);
    SimpleRead(std::string name, std::string seq, QualityValues qualities, SNR snr, Position qStart,
               Position qEnd, Frames pulseWidths);

    // general data
    std::string Name;
    std::string Sequence;
    QualityValues Qualities;
    SNR SignalToNoise;
    Position QueryStart;
    Position QueryEnd;
    boost::optional<Frames> PulseWidths;
};

class MappedSimpleRead : public SimpleRead
{
public:
    MappedSimpleRead(const SimpleRead& read, PacBio::BAM::Strand strand, Position templateStart,
                     Position templateEnd, PacBio::BAM::Cigar cigar, uint8_t mapQV);

    // mapping data
    PacBio::BAM::Strand Strand;
    Position TemplateStart;
    Position TemplateEnd;
    PacBio::BAM::Cigar Cigar;
    uint8_t MapQuality;
};

///
/// \brief
///
/// \param read
/// \param start
/// \param end
///
void ClipToQuery(SimpleRead& read, Position start, Position end);

///
/// \brief
///
/// \param read
/// \param start
/// \param end
///
void ClipToQuery(MappedSimpleRead& read, Position start, Position end);

///
/// \brief
///
/// \param read
/// \param start
/// \param end
/// \param exciseFlankingInserts
///
void ClipToReference(MappedSimpleRead& read, Position start, Position end,
                     bool exciseFlankingInserts);

}  // namespace BAM
}  // namespace PacBio

#endif  // SIMPLEREAD_H

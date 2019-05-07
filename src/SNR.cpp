// File Description
/// \file SNR.cpp
/// \brief Implements the SNR struct.
//
// Author: Lance Hepler, Derek Barnett

#include "PbbamInternalConfig.h"

#include "pbbam/SNR.h"

#include <algorithm>

#include <boost/algorithm/clamp.hpp>

namespace PacBio {
namespace BAM {

SNR::SNR(const double a, const double c, const double g, const double t) : A(a), C(c), G(g), T(t) {}

SNR::SNR(const std::vector<float>& snrs)
    : A(static_cast<double>(snrs[0]))
    , C(static_cast<double>(snrs[1]))
    , G(static_cast<double>(snrs[2]))
    , T(static_cast<double>(snrs[3]))
{
    assert(snrs.size() == 4);
}

SNR::SNR(const std::vector<double>& snrs) : A(snrs[0]), C(snrs[1]), G(snrs[2]), T(snrs[3])
{
    assert(snrs.size() == 4);
}

SNR::SNR(const double (&snrs)[4]) : A{snrs[0]}, C{snrs[1]}, G{snrs[2]}, T{snrs[3]} {}

SNR::operator std::vector<float>() const
{
    std::vector<float> snr = {static_cast<float>(A), static_cast<float>(C), static_cast<float>(G),
                              static_cast<float>(T)};
    return snr;
}

const double& SNR::operator[](const size_t i) const
{
    if (i == 0) return A;
    if (i == 1) return C;
    if (i == 2) return G;
    if (i == 3) return T;
    throw std::invalid_argument("SNR out of bounds!");
}

double& SNR::operator[](const size_t i)
{
    // casting away const when underlying object is non-const, is well-defined
    return const_cast<double&>(static_cast<const SNR&>(*this)[i]);
}

bool SNR::operator==(const SNR& other) const
{
    return std::tie(A, C, G, T) == std::tie(other.A, other.C, other.G, other.T);
}

bool SNR::operator!=(const SNR& other) const { return !(*this == other); }

double SNR::Minimum() const { return std::min(std::min(A, C), std::min(G, T)); }

SNR ClampSNR(const SNR& val, const SNR& lo, const SNR& hi)
{
    return SNR{
        boost::algorithm::clamp(val.A, lo.A, hi.A), boost::algorithm::clamp(val.C, lo.C, hi.C),
        boost::algorithm::clamp(val.G, lo.G, hi.G), boost::algorithm::clamp(val.T, lo.T, hi.T)};
}

}  // namespace BAM
}  // namespace PacBio

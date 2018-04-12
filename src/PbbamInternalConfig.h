// File Description
/// \file PbbamInternalConfig.h
/// \brief Defines internal macros for symbol visibility
//
// Author: Derek Barnett

#ifndef PBBAMINTERNALCONFIG_H
#define PBBAMINTERNALCONFIG_H

#if defined(WIN32)
#define PBBAM_EXPORT __declspec(dllexport)
#else
#define PBBAM_EXPORT __attribute__((visibility("default")))
#endif

#include "pbbam/Config.h"

#endif  // PBBAMINTERNALCONFIG_H

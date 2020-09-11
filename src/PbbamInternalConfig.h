#ifndef PBBAM_PBBAMINTERNALCONFIG_H
#define PBBAM_PBBAMINTERNALCONFIG_H

#if defined(WIN32)
#define PBBAM_EXPORT __declspec(dllexport)
#else
#define PBBAM_EXPORT __attribute__((visibility("default")))
#endif

#include <pbbam/Config.h>

#endif  // PBBAM_PBBAMINTERNALCONFIG_H

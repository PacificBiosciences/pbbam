/* VirtualPolymeraseReader.i */

%module PacBioBam

%{
#include <pbbam/virtual/VirtualPolymeraseReader.h>
#include <pbbam/virtual/ZmwReadStitcher.h>
using namespace PacBio;
using namespace PacBio::BAM;
typedef PacBio::BAM::ZmwReadStitcher VirtualPolymeraseReader;
%}

%include <pbbam/virtual/VirtualPolymeraseReader.h>
%include <pbbam/virtual/ZmwReadStitcher.h>
typedef PacBio::BAM::ZmwReadStitcher VirtualPolymeraseReader;

#ifdef SWIGPYTHON
%pythoncode %{

VirtualPolymeraseReader = ZmwReadStitcher

%}
#endif 
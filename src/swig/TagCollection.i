/* TagCollection.i */

%module PacBioBam

%{
#include <pbbam/TagCollection.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%template(TagCollectionType) std::map<std::string, PacBio::BAM::Tag>;

%include <pbbam/TagCollection.h>
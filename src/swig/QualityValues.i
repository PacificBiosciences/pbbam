/* QualityValues.i */
%module PacBioBam
%{
#include <pbbam/QualityValues.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%template(QualityValueList) std::vector<PacBio::BAM::QualityValue>;

%ignore PacBio::BAM::QualityValues::operator=;
%ignore PacBio::BAM::QualityValues::QualityValues(QualityValues&&); 

%include <pbbam/QualityValues.h>

/* CigarOperation.i */
%module PacBioBam
%{
#include <pbbam/CigarOperation.h>
using namespace PacBio;
using namespace PacBio::BAM;
%}

%ignore PacBio::BAM::CigarOperation::CigarOperation(CigarOperation&&);
%ignore PacBio::BAM::CigarOperation::operator=;

#ifdef SWIGR
%ignore PacBio::BAM::CigarOperation::CigarOperation(CigarOperationType, uint32_t);
#endif

%include <pbbam/CigarOperation.h>

// enums aren't always named consistently (at least between Mac/clang/swig & Linux/gcc/swig)
// so, keep this after the main %include so client sourcecan be consistent 
#ifdef SWIGPYTHON
%pythoncode %{
try:
	UNKNOWN_OP
	ALIGNMENT_MATCH
	INSERTION
	DELETION
	REFERENCE_SKIP
	SOFT_CLIP
	HARD_CLIP
	PADDING
	SEQUENCE_MATCH
	SEQUENCE_MISMATCH
except NameError:
	UNKNOWN_OP = CigarOperationType_UNKNOWN_OP
	ALIGNMENT_MATCH = CigarOperationType_ALIGNMENT_MATCH
	INSERTION = CigarOperationType_INSERTION
	DELETION = CigarOperationType_DELETION
	REFERENCE_SKIP = CigarOperationType_REFERENCE_SKIP
	SOFT_CLIP = CigarOperationType_SOFT_CLIP
	HARD_CLIP = CigarOperationType_HARD_CLIP
	PADDING = CigarOperationType_PADDING 
	SEQUENCE_MATCH = CigarOperationType_SEQUENCE_MATCH 
	SEQUENCE_MISMATCH = CigarOperationType_SEQUENCE_MISMATCH
%}
#endif


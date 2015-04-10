/* pbbam.i */
%module PacBioBam
%{
	
/*ifdef SWIGR
 define SWIG_SHARED_PTR_NAMESPACE boost
 define SWIG_SHARED_PTR_SUBNAMESPACE
endif*/
		
#include <pbbam/Config.h>
#include <string>
#include <vector>	
%}

#define SWIG_FILE_WITH_INIT
#define PBBAM_EXPORT

/********* SWIG includes ************/

%include "stdint.i"
%include "std_common.i"

#ifdef SWIGPYTHON
%include "std_shared_ptr.i"
#endif

#ifdef SWIGR
%include "boost_shared_ptr.i"
#endif

%include "std_container.i"
%include "std_map.i"
%include "std_pair.i"
%include "std_string.i"
%include "std_vector.i"

%template(StringList) std::vector<std::string>;

/********* PacBioBAM includes ************/

// Basic types
%include "Accuracy.i"
%include "CigarOperation.i"
%include "Interval.i"
%include "Orientation.i"
%include "Position.i"
%include "QualityValue.i"
%include "Strand.i"
%include "Tag.i"

// Basic type aggregates
%include "Cigar.i"
%include "GenomicInterval.i"
%include "QualityValues.i"
%include "TagCollection.i"

// keep this guy after the other basic types, hacky but works
%include "Frames.i" 

// Header API components
%include "ProgramInfo.i"
%include "ReadGroupInfo.i"
%include "SequenceInfo.i"
%include "BamHeader.i"

// SAM/BAM format
%include "BamFile.i"
%include "BamRecord.i"
%include "BamRecordImpl.i"
%include "BamTagCodec.i"
%include "BamWriter.i"
%include "SamTagCodec.i"

// Query/iterator API
%include "QueryBase.i"
%include "EntireFileQuery.i"
%include "GenomicIntervalQuery.i"

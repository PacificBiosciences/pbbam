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

#ifdef SWIGCSHARP
%rename(Equals)        *::operator==;
%rename(ToBool)        *::operator bool;
%rename(ToInt)         *::operator int;
%rename(ToUint8)       *::operator uint8_t;

%ignore *::operator !=;

// Iterator interfaces are not useful outside of C++
%ignore *::begin;
%ignore *::end;

%csmethodmodifiers *::ToString() const "public override";

#endif // SWIGCSHARP

/********* SWIG includes ************/

%include "stdint.i"
%include "std_common.i"

#ifdef SWIGR
%include "boost_shared_ptr.i"
#else
%include "std_shared_ptr.i"
#endif

%include "std_map.i"
%include "std_pair.i"
%include "std_string.i"
%include "std_vector.i"

 // TODO: can we call these vectors!?
%template(StringList) std::vector<std::string>;
%template(IntList)    std::vector<int>;
%template(UIntList)   std::vector<unsigned int>;
%template(FloatList)  std::vector<float>;
%template(ShortList)  std::vector<short>;
%template(CharList)   std::vector<char>;

// exception handling
%include "exception.i"
%exception {
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

/********* PacBioBAM includes ************/

#ifdef SWIGCSHARP
 // Renames to play nice with C#
 // (These are used in the dataset support code, where things like
 //  this happen in C++:
 //
 //    void Extensions(Extensions x) { ... }
 //
 //  and this poses problems for C#.  Renaming should be fine
 //  as it is doubtful we will refer to these classes by name anyway.)
 //
%rename(ExtensionsType)        PacBio::BAM::Extensions;
%rename(ExternalResourcesType) PacBio::BAM::ExternalResources;
%rename(FiltersType)           PacBio::BAM::Filters;
%rename(SubDataSetsType)       PacBio::BAM::SubDataSets;
%rename(ProvenanceType)        PacBio::BAM::Provenance;
%rename(PropertiesType)        PacBio::BAM::Properties;
%rename(FileIndicesType)       PacBio::BAM::FileIndices;
%rename(ParentToolType)        PacBio::BAM::ParentTool;
%rename(CigarType)             PacBio::BAM::Cigar;
#endif

// Basic types
%include "Accuracy.i"
%include "CigarOperation.i"
%include "Interval.i"
%include "LocalContextFlags.i"
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

// Header API 
%include "ProgramInfo.i"
%include "ReadGroupInfo.i"
%include "SequenceInfo.i"
%include "BamHeader.i"

// SAM/BAM format
%include "BamFile.i"
%include "BamRecordImpl.i"
%include "BamRecord.i"
%include "BamRecordBuilder.i"
%include "BamTagCodec.i"
%include "BamWriter.i"
%include "SamTagCodec.i"

// DataSet
%include "DataSetTypes.i"
%include "DataSet.i"

// Query/iterator API
%include "QueryBase.i"
%include "EntireFileQuery.i"
%include "GenomicIntervalQuery.i"
%include "ZmwQuery.i"
%include "ZmwGroupQuery.i"

// PBI
%include "PbiFile.i"
%include "PbiRawData.i"
%include "PbiIndex.i"

// FASTA
%include "IndexedFastaReader.i"

// VirtualPolymeraseBamRecord 
%include "VirtualRegion.i"
%include "VirtualPolymeraseBamRecord.i"
%include "VirtualPolymeraseReader.i"
%include "ZmwWhitelistVirtualReader.i"


# headers
set( PacBioBAM_H

    # API headers
    ${PacBioBAM_IncludeDir}/pbbam/Accuracy.h
    ${PacBioBAM_IncludeDir}/pbbam/BamFile.h
    ${PacBioBAM_IncludeDir}/pbbam/BamHeader.h
    ${PacBioBAM_IncludeDir}/pbbam/BamReader.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecordBuilder.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecordImpl.h
    ${PacBioBAM_IncludeDir}/pbbam/BamTagCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/BamWriter.h
    ${PacBioBAM_IncludeDir}/pbbam/Cigar.h
    ${PacBioBAM_IncludeDir}/pbbam/Config.h
    ${PacBioBAM_IncludeDir}/pbbam/EntireFileQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Frames.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicInterval.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicIntervalQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/IndexedFastaReader.h
    ${PacBioBAM_IncludeDir}/pbbam/Interval.h
    ${PacBioBAM_IncludeDir}/pbbam/Orientation.h
    ${PacBioBAM_IncludeDir}/pbbam/Position.h
    ${PacBioBAM_IncludeDir}/pbbam/ProgramInfo.h
    ${PacBioBAM_IncludeDir}/pbbam/QualityValue.h
    ${PacBioBAM_IncludeDir}/pbbam/QualityValues.h
    ${PacBioBAM_IncludeDir}/pbbam/QueryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/ReadGroupInfo.h
    ${PacBioBAM_IncludeDir}/pbbam/SamTagCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/SequenceInfo.h
    ${PacBioBAM_IncludeDir}/pbbam/Strand.h
    ${PacBioBAM_IncludeDir}/pbbam/Tag.h
    ${PacBioBAM_IncludeDir}/pbbam/TagCollection.h
#    ${PacBioBAM_IncludeDir}/pbbam/UnmappedReadsQuery.h

    # internal headers
    ${PacBioBAM_SourceDir}/AssertUtils.h
    ${PacBioBAM_SourceDir}/MemoryUtils.h
    ${PacBioBAM_SourceDir}/SequenceUtils.h
    ${PacBioBAM_SourceDir}/StringUtils.h
)

# sources
set( PacBioBAM_CPP

    ${PacBioBAM_SourceDir}/Accuracy.cpp
    ${PacBioBAM_SourceDir}/AssertUtils.cpp
    ${PacBioBAM_SourceDir}/BamFile.cpp
    ${PacBioBAM_SourceDir}/BamHeader.cpp
    ${PacBioBAM_SourceDir}/BamReader.cpp
    ${PacBioBAM_SourceDir}/BamRecord.cpp
    ${PacBioBAM_SourceDir}/BamRecordBuilder.cpp
    ${PacBioBAM_SourceDir}/BamRecordImpl.cpp
    ${PacBioBAM_SourceDir}/BamTagCodec.cpp
    ${PacBioBAM_SourceDir}/BamWriter.cpp
    ${PacBioBAM_SourceDir}/Cigar.cpp
    ${PacBioBAM_SourceDir}/CigarOperation.cpp
    ${PacBioBAM_SourceDir}/EntireFileQuery.cpp
    ${PacBioBAM_SourceDir}/Frames.cpp
    ${PacBioBAM_SourceDir}/GenomicInterval.cpp
    ${PacBioBAM_SourceDir}/GenomicIntervalQuery.cpp
    ${PacBioBAM_SourceDir}/IndexedFastaReader.cpp
    ${PacBioBAM_SourceDir}/MemoryUtils.cpp
    ${PacBioBAM_SourceDir}/ProgramInfo.cpp
    ${PacBioBAM_SourceDir}/QualityValue.cpp
    ${PacBioBAM_SourceDir}/QueryBase.cpp
    ${PacBioBAM_SourceDir}/ReadGroupInfo.cpp
    ${PacBioBAM_SourceDir}/SamTagCodec.cpp
    ${PacBioBAM_SourceDir}/SequenceInfo.cpp
    ${PacBioBAM_SourceDir}/Tag.cpp
    ${PacBioBAM_SourceDir}/TagCollection.cpp
#    ${PacBioBAM_SourceDir}/UnmappedReadsQuery.cpp
)

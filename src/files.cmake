
# headers
set( PacBioBAM_H

    # API headers
    ${PacBioBAM_IncludeDir}/pbbam/BamFile.h
    ${PacBioBAM_IncludeDir}/pbbam/BamReader.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/BamTagCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/BamWriter.h
    ${PacBioBAM_IncludeDir}/pbbam/Cigar.h
    ${PacBioBAM_IncludeDir}/pbbam/Config.h
    ${PacBioBAM_IncludeDir}/pbbam/DictionaryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/EntireFileQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicInterval.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicIntervalQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Interval.h
    ${PacBioBAM_IncludeDir}/pbbam/QueryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/SamHeader.h
    ${PacBioBAM_IncludeDir}/pbbam/SamHeaderCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/SamProgram.h
    ${PacBioBAM_IncludeDir}/pbbam/SamReadGroup.h
    ${PacBioBAM_IncludeDir}/pbbam/SamSequence.h
    ${PacBioBAM_IncludeDir}/pbbam/SamTagCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/Tag.h
    ${PacBioBAM_IncludeDir}/pbbam/TagCollection.h
    ${PacBioBAM_IncludeDir}/pbbam/UnmappedReadsQuery.h

    # internal headers
    ${PacBioBAM_SourceDir}/AssertUtils.h
    ${PacBioBAM_SourceDir}/MemoryUtils.h
)

# sources
set( PacBioBAM_CPP

    ${PacBioBAM_SourceDir}/AssertUtils.cpp
    ${PacBioBAM_SourceDir}/BamFile.cpp
    ${PacBioBAM_SourceDir}/BamReader.cpp
    ${PacBioBAM_SourceDir}/BamRecord.cpp
    ${PacBioBAM_SourceDir}/BamTagCodec.cpp
    ${PacBioBAM_SourceDir}/BamWriter.cpp
    ${PacBioBAM_SourceDir}/Cigar.cpp
    ${PacBioBAM_SourceDir}/EntireFileQuery.cpp
    ${PacBioBAM_SourceDir}/GenomicInterval.cpp
    ${PacBioBAM_SourceDir}/GenomicIntervalQuery.cpp
    ${PacBioBAM_SourceDir}/QueryBase.cpp
    ${PacBioBAM_SourceDir}/SamHeader.cpp
    ${PacBioBAM_SourceDir}/SamHeaderCodec.cpp
    ${PacBioBAM_SourceDir}/SamTagCodec.cpp
    ${PacBioBAM_SourceDir}/Tag.cpp
    ${PacBioBAM_SourceDir}/TagCollection.cpp
    ${PacBioBAM_SourceDir}/UnmappedReadsQuery.cpp

)

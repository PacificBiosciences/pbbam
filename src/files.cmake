
# headers
set( PacBioBAM_H

    # API headers
    ${PacBioBAM_IncludeDir}/pbbam/Accuracy.h
    ${PacBioBAM_IncludeDir}/pbbam/AlignmentPrinter.h
    ${PacBioBAM_IncludeDir}/pbbam/BamFile.h
    ${PacBioBAM_IncludeDir}/pbbam/BamHeader.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecordBuilder.h
    ${PacBioBAM_IncludeDir}/pbbam/BamRecordImpl.h
    ${PacBioBAM_IncludeDir}/pbbam/BamTagCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/BamWriter.h
    ${PacBioBAM_IncludeDir}/pbbam/Cigar.h
    ${PacBioBAM_IncludeDir}/pbbam/CigarOperation.h
    ${PacBioBAM_IncludeDir}/pbbam/Config.h
    ${PacBioBAM_IncludeDir}/pbbam/DataSet.h
    ${PacBioBAM_IncludeDir}/pbbam/DataSetTypes.h
    ${PacBioBAM_IncludeDir}/pbbam/EntireFileQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Frames.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicInterval.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicIntervalQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/GroupQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/GroupQueryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/IndexedFastaReader.h
    ${PacBioBAM_IncludeDir}/pbbam/Interval.h
    ${PacBioBAM_IncludeDir}/pbbam/LocalContextFlags.h
    ${PacBioBAM_IncludeDir}/pbbam/Orientation.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiBuilder.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiFile.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiIndex.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiRawData.h
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
    ${PacBioBAM_IncludeDir}/pbbam/ZmwGroupQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/ZmwQuery.h

    # internal headers
    ${PacBioBAM_IncludeDir}/pbbam/internal/BamRecordSort.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSet.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetBaseTypes.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetBaseTypes.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetElement.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetListElement.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetTypes.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/FilterEngine.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/IBamFileIterator.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/IMergeStrategy.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/MergeItem.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/MergeStrategy.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiIndex_p.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiIndex_p.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/QueryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/SequentialMergeStrategy.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/Tag.inl

    # virtual headers
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseBamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseReader.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegion.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegionType.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegionTypeMap.h

    ${PacBioBAM_SourceDir}/AssertUtils.h
    ${PacBioBAM_SourceDir}/DataSetIO.h
    ${PacBioBAM_SourceDir}/DataSetUtils.h
    ${PacBioBAM_SourceDir}/FileUtils.h
    ${PacBioBAM_SourceDir}/FofnReader.h
    ${PacBioBAM_SourceDir}/MemoryUtils.h
    ${PacBioBAM_SourceDir}/PbiIndexIO.h
    ${PacBioBAM_SourceDir}/SequenceUtils.h
    ${PacBioBAM_SourceDir}/StringUtils.h
    ${PacBioBAM_SourceDir}/TimeUtils.h
    ${PacBioBAM_SourceDir}/XmlReader.h
    ${PacBioBAM_SourceDir}/XmlWriter.h
    ${PacBioBAM_SourceDir}/pugixml/pugiconfig.hpp
    ${PacBioBAM_SourceDir}/pugixml/pugixml.hpp
)

# sources
set( PacBioBAM_CPP

    # main API headers
    ${PacBioBAM_SourceDir}/Accuracy.cpp
    ${PacBioBAM_SourceDir}/AlignmentPrinter.cpp
    ${PacBioBAM_SourceDir}/AssertUtils.cpp
    ${PacBioBAM_SourceDir}/BamFile.cpp
    ${PacBioBAM_SourceDir}/BamHeader.cpp
    ${PacBioBAM_SourceDir}/BamRecord.cpp
    ${PacBioBAM_SourceDir}/BamRecordBuilder.cpp
    ${PacBioBAM_SourceDir}/BamRecordImpl.cpp
    ${PacBioBAM_SourceDir}/BamTagCodec.cpp
    ${PacBioBAM_SourceDir}/BamWriter.cpp
    ${PacBioBAM_SourceDir}/Cigar.cpp
    ${PacBioBAM_SourceDir}/CigarOperation.cpp
    ${PacBioBAM_SourceDir}/Config.cpp
    ${PacBioBAM_SourceDir}/DataSet.cpp
    ${PacBioBAM_SourceDir}/DataSetBaseTypes.cpp
    ${PacBioBAM_SourceDir}/DataSetElement.cpp
    ${PacBioBAM_SourceDir}/DataSetIO.cpp
    ${PacBioBAM_SourceDir}/DataSetTypes.cpp
    ${PacBioBAM_SourceDir}/EntireFileQuery.cpp
    ${PacBioBAM_SourceDir}/FilterEngine.cpp
    ${PacBioBAM_SourceDir}/FofnReader.cpp
    ${PacBioBAM_SourceDir}/Frames.cpp
    ${PacBioBAM_SourceDir}/GenomicInterval.cpp
    ${PacBioBAM_SourceDir}/GenomicIntervalQuery.cpp
    ${PacBioBAM_SourceDir}/GroupQuery.cpp
    ${PacBioBAM_SourceDir}/IndexedFastaReader.cpp
    ${PacBioBAM_SourceDir}/MemoryUtils.cpp
    ${PacBioBAM_SourceDir}/PbiBuilder.cpp
    ${PacBioBAM_SourceDir}/PbiFile.cpp
    ${PacBioBAM_SourceDir}/PbiIndex.cpp
    ${PacBioBAM_SourceDir}/PbiIndexIO.cpp
    ${PacBioBAM_SourceDir}/PbiRawData.cpp
    ${PacBioBAM_SourceDir}/ProgramInfo.cpp
    ${PacBioBAM_SourceDir}/QualityValue.cpp
    ${PacBioBAM_SourceDir}/QueryBase.cpp
    ${PacBioBAM_SourceDir}/ReadGroupInfo.cpp
    ${PacBioBAM_SourceDir}/SamTagCodec.cpp
    ${PacBioBAM_SourceDir}/SequenceInfo.cpp
    ${PacBioBAM_SourceDir}/Tag.cpp
    ${PacBioBAM_SourceDir}/TagCollection.cpp
#    ${PacBioBAM_SourceDir}/UnmappedReadsQuery.cpp
    ${PacBioBAM_SourceDir}/XmlReader.cpp
    ${PacBioBAM_SourceDir}/XmlWriter.cpp
    ${PacBioBAM_SourceDir}/ZmwGroupQuery.cpp
    ${PacBioBAM_SourceDir}/ZmwQuery.cpp

    # virtual
    ${PacBioBAM_SourceDir}/VirtualPolymeraseBamRecord.cpp
    ${PacBioBAM_SourceDir}/VirtualPolymeraseReader.cpp
    ${PacBioBAM_SourceDir}/VirtualRegionTypeMap.cpp

    # XML I/O
    ${PacBioBAM_SourceDir}/pugixml/pugixml.cpp
)

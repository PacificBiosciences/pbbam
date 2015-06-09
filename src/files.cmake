
# headers
set( PacBioBAM_H

    # API headers
    ${PacBioBAM_IncludeDir}/pbbam/Accuracy.h
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

    # dataset headers
    ${PacBioBAM_IncludeDir}/pbbam/dataset/AlignmentSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/AlignmentSetMetadata.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/BarcodeSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/BarcodeSetMetadata.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/CcsReadSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/CcsReadSetMetadata.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/ContigSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/ContigSetMetadata.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/DataSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/DataSetBase.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/DataSetMetadata.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/DataSetMetadataBase.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/ExternalDataReferences.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/Filters.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/ReferenceSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/ReferenceSetMetadata.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/SubDataSets.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/SubreadSet.h
    ${PacBioBAM_IncludeDir}/pbbam/dataset/SubreadSetMetadata.h

    # internal headers
    ${PacBioBAM_IncludeDir}/pbbam/internal/BamRecordSort.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetElement.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetListElement.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/FilterEngine.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/IBamFileIterator.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/IMergeStrategy.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/MergeStrategy.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/QueryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/MergeItem.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/SequentialMergeStrategy.h

    # virtual headers
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseBamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseReader.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegion.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegionType.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegionTypeMap.h

    ${PacBioBAM_SourceDir}/AssertUtils.h
    ${PacBioBAM_SourceDir}/DataSetIO.h
    ${PacBioBAM_SourceDir}/FileUtils.h
    ${PacBioBAM_SourceDir}/FofnReader.h
    ${PacBioBAM_SourceDir}/MemoryUtils.h
    ${PacBioBAM_SourceDir}/PbiIndex_p.h
    ${PacBioBAM_SourceDir}/PbiIndexIO.h
    ${PacBioBAM_SourceDir}/SequenceUtils.h
    ${PacBioBAM_SourceDir}/StringUtils.h
    ${PacBioBAM_SourceDir}/XmlReader.h
    ${PacBioBAM_SourceDir}/XmlWriter.h
    ${PacBioBAM_SourceDir}/pugixml/pugiconfig.hpp
    ${PacBioBAM_SourceDir}/pugixml/pugixml.hpp
)

# sources
set( PacBioBAM_CPP

    # main API headers
    ${PacBioBAM_SourceDir}/Accuracy.cpp
    ${PacBioBAM_SourceDir}/AlignmentSet.cpp
    ${PacBioBAM_SourceDir}/AlignmentSetMetadata.cpp
    ${PacBioBAM_SourceDir}/AssertUtils.cpp
    ${PacBioBAM_SourceDir}/BamFile.cpp
    ${PacBioBAM_SourceDir}/BamHeader.cpp
    ${PacBioBAM_SourceDir}/BamRecord.cpp
    ${PacBioBAM_SourceDir}/BamRecordBuilder.cpp
    ${PacBioBAM_SourceDir}/BamRecordImpl.cpp
    ${PacBioBAM_SourceDir}/BamTagCodec.cpp
    ${PacBioBAM_SourceDir}/BamWriter.cpp
    ${PacBioBAM_SourceDir}/BarcodeSet.cpp
    ${PacBioBAM_SourceDir}/BarcodeSetMetadata.cpp
    ${PacBioBAM_SourceDir}/CcsReadSet.cpp
    ${PacBioBAM_SourceDir}/CcsReadSetMetadata.cpp
    ${PacBioBAM_SourceDir}/Cigar.cpp
    ${PacBioBAM_SourceDir}/CigarOperation.cpp
    ${PacBioBAM_SourceDir}/Config.cpp
    ${PacBioBAM_SourceDir}/ContigSet.cpp
    ${PacBioBAM_SourceDir}/ContigSetMetadata.cpp
    ${PacBioBAM_SourceDir}/DataSet.cpp
    ${PacBioBAM_SourceDir}/DataSetBase.cpp
    ${PacBioBAM_SourceDir}/DataSetElement.cpp
    ${PacBioBAM_SourceDir}/DataSetIO.cpp
    ${PacBioBAM_SourceDir}/DataSetMetadata.cpp
    ${PacBioBAM_SourceDir}/DataSetMetadataBase.cpp
    ${PacBioBAM_SourceDir}/EntireFileQuery.cpp
    ${PacBioBAM_SourceDir}/ExternalDataReferences.cpp
    ${PacBioBAM_SourceDir}/FilterEngine.cpp
    ${PacBioBAM_SourceDir}/Filters.cpp
    ${PacBioBAM_SourceDir}/FofnReader.cpp
    ${PacBioBAM_SourceDir}/Frames.cpp
    ${PacBioBAM_SourceDir}/GenomicInterval.cpp
    ${PacBioBAM_SourceDir}/GenomicIntervalQuery.cpp
    ${PacBioBAM_SourceDir}/GroupQuery.cpp
    ${PacBioBAM_SourceDir}/IndexedFastaReader.cpp
    ${PacBioBAM_SourceDir}/MemoryUtils.cpp
    ${PacBioBAM_SourceDir}/PbiFile.cpp
    ${PacBioBAM_SourceDir}/PbiIndex.cpp
    ${PacBioBAM_SourceDir}/PbiIndexIO.cpp
    ${PacBioBAM_SourceDir}/PbiRawData.cpp
    ${PacBioBAM_SourceDir}/ProgramInfo.cpp
    ${PacBioBAM_SourceDir}/QualityValue.cpp
    ${PacBioBAM_SourceDir}/QueryBase.cpp
    ${PacBioBAM_SourceDir}/ReadGroupInfo.cpp
    ${PacBioBAM_SourceDir}/ReferenceSet.cpp
    ${PacBioBAM_SourceDir}/ReferenceSetMetadata.cpp
    ${PacBioBAM_SourceDir}/SamTagCodec.cpp
    ${PacBioBAM_SourceDir}/SequenceInfo.cpp
    ${PacBioBAM_SourceDir}/SubDataSets.cpp
    ${PacBioBAM_SourceDir}/SubreadSet.cpp
    ${PacBioBAM_SourceDir}/SubreadSetMetadata.cpp
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


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
    ${PacBioBAM_IncludeDir}/pbbam/BaiIndexedBamReader.h
    ${PacBioBAM_IncludeDir}/pbbam/BamReader.h
    ${PacBioBAM_IncludeDir}/pbbam/CompositeBamReader.h
    ${PacBioBAM_IncludeDir}/pbbam/BamWriter.h
    ${PacBioBAM_IncludeDir}/pbbam/BarcodeQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Cigar.h
    ${PacBioBAM_IncludeDir}/pbbam/CigarOperation.h
    ${PacBioBAM_IncludeDir}/pbbam/Compare.h
    ${PacBioBAM_IncludeDir}/pbbam/Config.h
    ${PacBioBAM_IncludeDir}/pbbam/DataSet.h
    ${PacBioBAM_IncludeDir}/pbbam/DataSetTypes.h
    ${PacBioBAM_IncludeDir}/pbbam/DataSetXsd.h
    ${PacBioBAM_IncludeDir}/pbbam/EntireFileQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Frames.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicInterval.h
    ${PacBioBAM_IncludeDir}/pbbam/GenomicIntervalQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/IndexedFastaReader.h
    ${PacBioBAM_IncludeDir}/pbbam/Interval.h
    ${PacBioBAM_IncludeDir}/pbbam/LocalContextFlags.h
    ${PacBioBAM_IncludeDir}/pbbam/MD5.h
    ${PacBioBAM_IncludeDir}/pbbam/Orientation.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiBasicTypes.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiBuilder.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiFile.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiFilter.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiFilterQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiFilterTypes.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiIndex.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiIndexedBamReader.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiLookupData.h
    ${PacBioBAM_IncludeDir}/pbbam/PbiRawData.h
    ${PacBioBAM_IncludeDir}/pbbam/Position.h
    ${PacBioBAM_IncludeDir}/pbbam/ProgramInfo.h
    ${PacBioBAM_IncludeDir}/pbbam/QNameQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/QualityValue.h
    ${PacBioBAM_IncludeDir}/pbbam/QualityValues.h
    ${PacBioBAM_IncludeDir}/pbbam/ReadAccuracyQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/ReadGroupInfo.h
    ${PacBioBAM_IncludeDir}/pbbam/SamTagCodec.h
    ${PacBioBAM_IncludeDir}/pbbam/SequenceInfo.h
    ${PacBioBAM_IncludeDir}/pbbam/Strand.h  
    ${PacBioBAM_IncludeDir}/pbbam/SubreadLengthQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Tag.h
    ${PacBioBAM_IncludeDir}/pbbam/TagCollection.h
#    ${PacBioBAM_IncludeDir}/pbbam/UnmappedReadsQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/Validator.h
    ${PacBioBAM_IncludeDir}/pbbam/ZmwGroupQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/ZmwQuery.h
    ${PacBioBAM_IncludeDir}/pbbam/ZmwType.h
    ${PacBioBAM_IncludeDir}/pbbam/ZmwTypeMap.h

    # exception headers
    ${PacBioBAM_IncludeDir}/pbbam/exception/InvalidSequencingChemistryException.h
    ${PacBioBAM_IncludeDir}/pbbam/exception/ValidationException.h

    # API-internal headers & inline files
    ${PacBioBAM_IncludeDir}/pbbam/internal/Accuracy.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/BamHeader.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/BamRecord.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/BamRecordBuilder.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/BamRecordImpl.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/Cigar.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/CigarOperation.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/Compare.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/CompositeBamReader.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSet.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetBaseTypes.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetBaseTypes.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetElement.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetElement.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetListElement.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetListElement.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/DataSetTypes.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/Frames.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/GenomicInterval.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/Interval.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiBasicTypes.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiFilter.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiFilterTypes.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiIndex.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiLookupData.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/PbiRawData.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/ProgramInfo.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/QualityValue.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/QualityValues.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/QueryBase.h
    ${PacBioBAM_IncludeDir}/pbbam/internal/QueryBase.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/ReadGroupInfo.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/SequenceInfo.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/Tag.inl
    ${PacBioBAM_IncludeDir}/pbbam/internal/Validator.inl

    # virtual headers
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseBamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseCompositeReader.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualPolymeraseReader.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegion.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegionType.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualRegionTypeMap.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/VirtualZmwBamRecord.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/WhitelistedZmwReadStitcher.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/ZmwReadStitcher.h
    ${PacBioBAM_IncludeDir}/pbbam/virtual/ZmwWhitelistVirtualReader.h

    # library-internal headers
    ${PacBioBAM_SourceDir}/AssertUtils.h
    ${PacBioBAM_SourceDir}/ChemistryTable.h
    ${PacBioBAM_SourceDir}/DataSetIO.h
    ${PacBioBAM_SourceDir}/DataSetUtils.h
    ${PacBioBAM_SourceDir}/FileProducer.h
    ${PacBioBAM_SourceDir}/FileUtils.h
    ${PacBioBAM_SourceDir}/FofnReader.h
    ${PacBioBAM_SourceDir}/MemoryUtils.h
    ${PacBioBAM_SourceDir}/PbiIndexIO.h
    ${PacBioBAM_SourceDir}/SequenceUtils.h
    ${PacBioBAM_SourceDir}/StringUtils.h
    ${PacBioBAM_SourceDir}/TimeUtils.h
    ${PacBioBAM_SourceDir}/ValidationErrors.h
    ${PacBioBAM_SourceDir}/Version.h
    ${PacBioBAM_SourceDir}/VirtualZmwCompositeReader.h
    ${PacBioBAM_SourceDir}/VirtualZmwReader.h
    ${PacBioBAM_SourceDir}/XmlReader.h
    ${PacBioBAM_SourceDir}/XmlWriter.h
    ${PacBioBAM_SourceDir}/pugixml/pugiconfig.hpp
    ${PacBioBAM_SourceDir}/pugixml/pugixml.hpp
)

# sources
set( PacBioBAM_CPP

    ${PacBioBAM_SourceDir}/Accuracy.cpp
    ${PacBioBAM_SourceDir}/AlignmentPrinter.cpp
    ${PacBioBAM_SourceDir}/AssertUtils.cpp
    ${PacBioBAM_SourceDir}/BaiIndexedBamReader.cpp
    ${PacBioBAM_SourceDir}/BamFile.cpp
    ${PacBioBAM_SourceDir}/BamHeader.cpp
    ${PacBioBAM_SourceDir}/BamReader.cpp
    ${PacBioBAM_SourceDir}/BamRecord.cpp
    ${PacBioBAM_SourceDir}/BamRecordBuilder.cpp
    ${PacBioBAM_SourceDir}/BamRecordImpl.cpp
    ${PacBioBAM_SourceDir}/BamTagCodec.cpp
    ${PacBioBAM_SourceDir}/BamWriter.cpp
    ${PacBioBAM_SourceDir}/BarcodeQuery.cpp
    ${PacBioBAM_SourceDir}/ChemistryTable.cpp
    ${PacBioBAM_SourceDir}/Cigar.cpp
    ${PacBioBAM_SourceDir}/CigarOperation.cpp
    ${PacBioBAM_SourceDir}/Compare.cpp
    ${PacBioBAM_SourceDir}/Config.cpp
    ${PacBioBAM_SourceDir}/DataSet.cpp
    ${PacBioBAM_SourceDir}/DataSetBaseTypes.cpp
    ${PacBioBAM_SourceDir}/DataSetElement.cpp
    ${PacBioBAM_SourceDir}/DataSetIO.cpp
    ${PacBioBAM_SourceDir}/DataSetTypes.cpp
    ${PacBioBAM_SourceDir}/DataSetXsd.cpp
    ${PacBioBAM_SourceDir}/EntireFileQuery.cpp
    ${PacBioBAM_SourceDir}/FileProducer.cpp
    ${PacBioBAM_SourceDir}/FileUtils.cpp
    ${PacBioBAM_SourceDir}/FofnReader.cpp
    ${PacBioBAM_SourceDir}/Frames.cpp
    ${PacBioBAM_SourceDir}/GenomicInterval.cpp
    ${PacBioBAM_SourceDir}/GenomicIntervalQuery.cpp
    ${PacBioBAM_SourceDir}/IndexedFastaReader.cpp
    ${PacBioBAM_SourceDir}/MD5.cpp
    ${PacBioBAM_SourceDir}/MemoryUtils.cpp
    ${PacBioBAM_SourceDir}/PbiBuilder.cpp
    ${PacBioBAM_SourceDir}/PbiFile.cpp
    ${PacBioBAM_SourceDir}/PbiFilter.cpp
    ${PacBioBAM_SourceDir}/PbiFilterQuery.cpp
    ${PacBioBAM_SourceDir}/PbiFilterTypes.cpp
    ${PacBioBAM_SourceDir}/PbiIndex.cpp
    ${PacBioBAM_SourceDir}/PbiIndexedBamReader.cpp
    ${PacBioBAM_SourceDir}/PbiIndexIO.cpp
    ${PacBioBAM_SourceDir}/PbiRawData.cpp
    ${PacBioBAM_SourceDir}/ProgramInfo.cpp
    ${PacBioBAM_SourceDir}/QNameQuery.cpp
    ${PacBioBAM_SourceDir}/QualityValue.cpp
    ${PacBioBAM_SourceDir}/ReadAccuracyQuery.cpp
    ${PacBioBAM_SourceDir}/ReadGroupInfo.cpp
    ${PacBioBAM_SourceDir}/SamTagCodec.cpp
    ${PacBioBAM_SourceDir}/SequenceInfo.cpp
    ${PacBioBAM_SourceDir}/SubreadLengthQuery.cpp
    ${PacBioBAM_SourceDir}/Tag.cpp
    ${PacBioBAM_SourceDir}/TagCollection.cpp
#    ${PacBioBAM_SourceDir}/UnmappedReadsQuery.cpp
    ${PacBioBAM_SourceDir}/Validator.cpp
    ${PacBioBAM_SourceDir}/ValidationErrors.cpp
    ${PacBioBAM_SourceDir}/ValidationException.cpp
    ${PacBioBAM_SourceDir}/Version.cpp
    ${PacBioBAM_SourceDir}/VirtualZmwBamRecord.cpp
    ${PacBioBAM_SourceDir}/VirtualZmwCompositeReader.cpp
    ${PacBioBAM_SourceDir}/VirtualZmwReader.cpp
    ${PacBioBAM_SourceDir}/VirtualRegionTypeMap.cpp
    ${PacBioBAM_SourceDir}/XmlReader.cpp
    ${PacBioBAM_SourceDir}/XmlWriter.cpp
    ${PacBioBAM_SourceDir}/WhitelistedZmwReadStitcher.cpp
    ${PacBioBAM_SourceDir}/ZmwGroupQuery.cpp
    ${PacBioBAM_SourceDir}/ZmwReadStitcher.cpp
    ${PacBioBAM_SourceDir}/ZmwQuery.cpp
    ${PacBioBAM_SourceDir}/ZmwTypeMap.cpp

    # XML I/O
    ${PacBioBAM_SourceDir}/pugixml/pugixml.cpp
)

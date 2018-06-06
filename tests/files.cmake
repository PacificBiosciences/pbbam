# test case headers
set( PacBioBAMTest_H

)

# test case sources
set( PacBioBAMTest_CPP

    ${PacBioBAM_TestsDir}/src/test_Accuracy.cpp
    ${PacBioBAM_TestsDir}/src/test_AlignmentPrinter.cpp
    ${PacBioBAM_TestsDir}/src/test_BamFile.cpp
    ${PacBioBAM_TestsDir}/src/test_BamHeader.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecord.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordBuilder.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordClipping.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordImplCore.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordImplTags.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordImplVariableData.cpp
    ${PacBioBAM_TestsDir}/src/test_BamRecordMapping.cpp
    ${PacBioBAM_TestsDir}/src/test_BamWriter.cpp
    ${PacBioBAM_TestsDir}/src/test_BarcodeQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_Cigar.cpp
    ${PacBioBAM_TestsDir}/src/test_Compare.cpp
    ${PacBioBAM_TestsDir}/src/test_DataSetCore.cpp
    ${PacBioBAM_TestsDir}/src/test_DataSetIO.cpp
    ${PacBioBAM_TestsDir}/src/test_DataSetQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_DataSetXsd.cpp
    ${PacBioBAM_TestsDir}/src/test_EndToEnd.cpp
    ${PacBioBAM_TestsDir}/src/test_EntireFileQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_Fasta.cpp
    ${PacBioBAM_TestsDir}/src/test_Fastq.cpp
    ${PacBioBAM_TestsDir}/src/test_FileUtils.cpp
    ${PacBioBAM_TestsDir}/src/test_Frames.cpp
    ${PacBioBAM_TestsDir}/src/test_GenomicIntervalQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_IndexedBamWriter.cpp
    ${PacBioBAM_TestsDir}/src/test_IndexedFastaReader.cpp
    ${PacBioBAM_TestsDir}/src/test_Intervals.cpp
    ${PacBioBAM_TestsDir}/src/test_LongCigar.cpp
    ${PacBioBAM_TestsDir}/src/test_PacBioIndex.cpp
    ${PacBioBAM_TestsDir}/src/test_PbiFilter.cpp
    ${PacBioBAM_TestsDir}/src/test_PbiFilterQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_QNameQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_QualityValues.cpp
    ${PacBioBAM_TestsDir}/src/test_Pulse2BaseCache.cpp
    ${PacBioBAM_TestsDir}/src/test_ReadAccuracyQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_ReadGroupInfo.cpp
    ${PacBioBAM_TestsDir}/src/test_SamWriter.cpp
    ${PacBioBAM_TestsDir}/src/test_SequenceUtils.cpp
    ${PacBioBAM_TestsDir}/src/test_StringUtils.cpp
    ${PacBioBAM_TestsDir}/src/test_SubreadLengthQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_Tags.cpp
    ${PacBioBAM_TestsDir}/src/test_TimeUtils.cpp
    ${PacBioBAM_TestsDir}/src/test_Validator.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfFile.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfFormat.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfHeader.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfReader.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfSort.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfQuery.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfVariant.cpp
    ${PacBioBAM_TestsDir}/src/test_VcfWriter.cpp
    ${PacBioBAM_TestsDir}/src/test_Version.cpp
    ${PacBioBAM_TestsDir}/src/test_WhitelistedZmwReadStitcher.cpp
    ${PacBioBAM_TestsDir}/src/test_ZmwReadStitcher.cpp
    ${PacBioBAM_TestsDir}/src/test_ZmwQuery.cpp
)

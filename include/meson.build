###########
# headers #
###########

if not meson.is_subproject()
  install_headers(
    files([
      'pbbam/AlignmentPrinter.h',
      'pbbam/BaiIndexCache.h',
      'pbbam/BaiIndexedBamReader.h',
      'pbbam/BamFile.h',
      'pbbam/BamFileMerger.h',
      'pbbam/BamHeader.h',
      'pbbam/BamReader.h',
      'pbbam/BamRecord.h',
      'pbbam/BamRecordImpl.h',
      'pbbam/BamRecordTag.h',
      'pbbam/BamRecordView.h',
      'pbbam/BamTagCodec.h',
      'pbbam/BamWriter.h',
      'pbbam/BgzipFastaWriter.h',
      'pbbam/BgzipFastqWriter.h',
      'pbbam/BgzipWriter.h',
      'pbbam/ClipType.h',
      'pbbam/CollectionMetadata.h',
      'pbbam/Compare.h',
      'pbbam/CompositeBamReader.h',
      'pbbam/CompositeFastaReader.h',
      'pbbam/Config.h',
      'pbbam/DataSet.h',
      'pbbam/DataSetTypes.h',
      'pbbam/DataSetXsd.h',
      'pbbam/Deleters.h',
      'pbbam/EntireFileQuery.h',
      'pbbam/FaiIndex.h',
      'pbbam/FastaCache.h',
      'pbbam/FastaReader.h',
      'pbbam/FastaSequence.h',
      'pbbam/FastaSequenceQuery.h',
      'pbbam/FastaWriter.h',
      'pbbam/FastqReader.h',
      'pbbam/FastqSequence.h',
      'pbbam/FastqWriter.h',
      'pbbam/FormatUtils.h',
      'pbbam/FrameEncodingType.h',
      'pbbam/GenomicIntervalQuery.h',
      'pbbam/IFastaWriter.h',
      'pbbam/IFastqWriter.h',
      'pbbam/IndexedBamWriter.h',
      'pbbam/IndexedFastaReader.h',
      'pbbam/IndexedFastqReader.h',
      'pbbam/IRecordWriter.h',
      'pbbam/LibraryInfo.h',
      'pbbam/MD5.h',
      'pbbam/MoveAppend.h',
      'pbbam/PbbamVersion.h',
      'pbbam/PbiBasicTypes.h',
      'pbbam/PbiBuilder.h',
      'pbbam/PbiFile.h',
      'pbbam/PbiFilter.h',
      'pbbam/PbiFilterQuery.h',
      'pbbam/PbiFilterTypes.h',
      'pbbam/PbiIndexedBamReader.h',
      'pbbam/PbiRawData.h',
      'pbbam/ProgramInfo.h',
      'pbbam/PulseBehavior.h',
      'pbbam/PulseExclusionReason.h',
      'pbbam/ReadGroupInfo.h',
      'pbbam/RecordType.h',
      'pbbam/RunMetadata.h',
      'pbbam/SamReader.h',
      'pbbam/SamTagCodec.h',
      'pbbam/SamWriter.h',
      'pbbam/SequenceInfo.h',
      'pbbam/StringUtilities.h',
      'pbbam/Tag.h',
      'pbbam/TagCollection.h',
      'pbbam/TextFileReader.h',
      'pbbam/TextFileWriter.h',
      'pbbam/Validator.h',
      'pbbam/ZmwGroupQuery.h',
      'pbbam/ZmwQuery.h',
      'pbbam/ZmwType.h',
      'pbbam/ZmwTypeMap.h']),
    subdir : 'pbbam')

  install_headers(
    files([
      'pbbam/bed/BedReader.h',
      'pbbam/bed/BedWriter.h']),
    subdir : 'pbbam/bed')

  install_headers(
    files([
      'pbbam/csv/CsvReader.h',
      'pbbam/csv/CsvTypes.h',
      'pbbam/csv/CsvWriter.h']),
    subdir : 'pbbam/csv')

  install_headers(
    files([
      'pbbam/exception/BundleChemistryMappingException.h',
      'pbbam/exception/InvalidSequencingChemistryException.h',
      'pbbam/exception/ValidationException.h']),
    subdir : 'pbbam/exception')

  install_headers(
    files([
      'pbbam/internal/Compare.inl',
      'pbbam/internal/CompositeBamReader.inl',
      'pbbam/internal/DataSetBaseTypes.h',
      'pbbam/internal/DataSetElement.h',
      'pbbam/internal/DataSetElement.inl',
      'pbbam/internal/PbiBasicTypes.inl',
      'pbbam/internal/PbiFilter.inl',
      'pbbam/internal/PbiFilterTypes.inl',
      'pbbam/internal/QueryBase.h',
      'pbbam/internal/QueryBase.inl']),
    subdir : 'pbbam/internal')

  install_headers(
    files([
      'pbbam/vcf/VcfVariant.h',
      'pbbam/vcf/VcfFile.h',
      'pbbam/vcf/VcfFormat.h',
      'pbbam/vcf/VcfHeader.h',
      'pbbam/vcf/VcfHeaderTypes.h',
      'pbbam/vcf/VcfReader.h',
      'pbbam/vcf/VcfSort.h',
      'pbbam/vcf/VcfQuery.h',
      'pbbam/vcf/VcfWriter.h']),
    subdir : 'pbbam/vcf')

  install_headers(
    files([
      'pbbam/virtual/VirtualPolymeraseBamRecord.h',
      'pbbam/virtual/VirtualPolymeraseCompositeReader.h',
      'pbbam/virtual/VirtualPolymeraseReader.h',
      'pbbam/virtual/VirtualRegion.h',
      'pbbam/virtual/VirtualRegionType.h',
      'pbbam/virtual/VirtualRegionTypeMap.h',
      'pbbam/virtual/VirtualZmwBamRecord.h',
      'pbbam/virtual/WhitelistedZmwReadStitcher.h',
      'pbbam/virtual/ZmwReadStitcher.h',
      'pbbam/virtual/ZmwWhitelistVirtualReader.h']),
    subdir : 'pbbam/virtual')
endif

pbbam_include_directories = include_directories('.')

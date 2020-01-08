###########
# sources #
###########

pbbam_cpp_sources = files([
  'AlignmentPrinter.cpp',
  'BaiIndexCache.cpp',
  'BaiIndexedBamReader.cpp',
  'BamFile.cpp',
  'BamFileMerger.cpp',
  'BamHeader.cpp',
  'BamReader.cpp',
  'BamRecord.cpp',
  'BamRecordBuilder.cpp',
  'BamRecordImpl.cpp',
  'BamRecordTags.cpp',
  'BamRecordView.cpp',
  'BamTagCodec.cpp',
  'BamWriter.cpp',
  'BgzipFastaWriter.cpp',
  'BgzipFastqWriter.cpp',
  'BgzipWriter.cpp',
  'ChemistryTable.cpp',
  'Compare.cpp',
  'CompositeBamReader.cpp',
  'CompositeFastaReader.cpp',
  'Config.cpp',
  'DataSet.cpp',
  'DataSetBaseTypes.cpp',
  'DataSetElement.cpp',
  'DataSetIO.cpp',
  'DataSetTypes.cpp',
  'DataSetXsd.cpp',
  'EntireFileQuery.cpp',
  'FaiIndex.cpp',
  'FaiZmwChunker.cpp',
  'FastaCache.cpp',
  'FastaReader.cpp',
  'FastaSequence.cpp',
  'FastaSequenceQuery.cpp',
  'FastaWriter.cpp',
  'FastqReader.cpp',
  'FastqSequence.cpp',
  'FastqWriter.cpp',
  'FileProducer.cpp',
  'FileUtils.cpp',
  'FofnReader.cpp',
  'FormatUtils.cpp',
  'GenomicIntervalQuery.cpp',
  'IFastaWriter.cpp',
  'IFastqWriter.cpp',
  'IndexedBamWriter.cpp',
  'IndexedFastaReader.cpp',
  'IndexedFastqBgzfReader.cpp',
  'IndexedFastqReader.cpp',
  'IndexedFastqTextReader.cpp',
  'IndexedFastqReaderImpl.cpp',
  'IRecordWriter.cpp',
  'KSeqReader.cpp',
  'MD5.cpp',
  'MemoryUtils.cpp',
  'PbiBuilder.cpp',
  'PbiFile.cpp',
  'PbiFilter.cpp',
  'PbiFilterQuery.cpp',
  'PbiFilterTypes.cpp',
  'PbiIndexedBamReader.cpp',
  'PbiIndexIO.cpp',
  'PbiRawData.cpp',
  'ProgramInfo.cpp',
  'ReadGroupInfo.cpp',
  'RecordType.cpp',
  'RunMetadata.cpp',
  'RunMetadataParser.cpp',
  'SamTagCodec.cpp',
  'SamWriter.cpp',
  'SequenceInfo.cpp',
  'StringUtilities.cpp',
  'Tag.cpp',
  'TagCollection.cpp',
  'TextFileReader.cpp',
  'TextFileWriter.cpp',
  'Validator.cpp',
  'ValidationErrors.cpp',
  'ValidationException.cpp',
  'Version.cpp',
  'VirtualZmwBamRecord.cpp',
  'VirtualZmwCompositeReader.cpp',
  'VirtualZmwReader.cpp',
  'VirtualRegion.cpp',
  'VirtualRegionTypeMap.cpp',
  'XmlReader.cpp',
  'XmlWriter.cpp',
  'WhitelistedZmwReadStitcher.cpp',
  'ZmwChunkedFastaReader.cpp',
  'ZmwChunkedFastqReader.cpp',
  'ZmwChunkedFastxBgzfReader.cpp',
  'ZmwChunkedFastxReaderImpl.cpp',
  'ZmwChunkedFastxTextReader.cpp',
  'ZmwGroupQuery.cpp',
  'ZmwReadStitcher.cpp',
  'ZmwQuery.cpp',
  'ZmwTypeMap.cpp',

  # bed
  'bed/BedReader.cpp',
  'bed/BedWriter.cpp',

  # ccs
  'ccs/CCSPbiBuilder.cpp',
  'ccs/CCSRecord.cpp',
  'ccs/CCSRecordFormat.cpp',
  'ccs/CCSRecordReader.cpp',
  'ccs/CCSRecordWriter.cpp',

  # vcf
  'vcf/VcfFile.cpp',
  'vcf/VcfFormat.cpp',
  'vcf/VcfHeader.cpp',
  'vcf/VcfHeaderTypes.cpp',
  'vcf/VcfQuery.cpp',
  'vcf/VcfReader.cpp',
  'vcf/VcfSort.cpp',
  'vcf/VcfVariant.cpp',
  'vcf/VcfWriter.cpp',

  # XML I/O
  'pugixml/pugixml.cpp'
])

pbbam_extra_flags = []
if get_option('auto-validate')
  pbbam_extra_flags += '-DPBBAM_AUTOVALIDATE=1'
endif

# install library if
# - either running as a proper project
# - or using shared libraries
pbbam_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

pbbam_lib = library(
  'pbbam',
  pbbam_cpp_sources,
  # use boost SONAME practice:
  #   cause ld startup issues before
  #   you even have the chance of running
  #   into ABI issues.
  soversion : meson.project_version(),
  version : meson.project_version(),
  install : pbbam_lib_install,
  dependencies : [pbbam_thread_dep, pbbam_boost_dep, pbbam_zlib_dep, pbbam_htslib_dep, pbbam_pbcopper_dep],
  include_directories : pbbam_include_directories,
  cpp_args : [pbbam_extra_flags, pbbam_warning_flags, pbbam_macros])

################
# version info #
################

pbbam_ver_arr = meson.project_version().split('.')
pbbam_major_version = pbbam_ver_arr[0]
pbbam_minor_version = pbbam_ver_arr[1]
pbbam_patch_version = pbbam_ver_arr[2]

pbbam_config = configuration_data()
pbbam_config.set('PBBAM_VERSION', meson.project_version())
pbbam_config.set('PBBAM_VERSION_MAJOR', pbbam_major_version)
pbbam_config.set('PBBAM_VERSION_MINOR', pbbam_minor_version)
pbbam_config.set('PBBAM_VERSION_PATCH', pbbam_patch_version)

pbbam_version_sources = [

  # old <-- to be removed
  vcs_tag(
    input : 'PbbamGitHash.cpp.in',
    output : 'PbbamGitHash.cpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PBBAM_GIT_SHA1@'),
  configure_file(
    input : 'PbbamVersion.cpp.in',
    output : 'PbbamVersion.cpp',
    configuration : pbbam_config),

  # new
  vcs_tag(
    input : 'LibraryGitHash.h.in',
    output : 'LibraryGitHash.h',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PBBAM_GIT_SHA1@'),
  configure_file(
    input : 'LibraryVersion.h.in',
    output : 'LibraryVersion.h',
    configuration : pbbam_config)
]

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
  'BamRecordImpl.cpp',
  'BamRecordTags.cpp',
  'BamRecordView.cpp',
  'BamTagCodec.cpp',
  'BamWriter.cpp',
  'BgzipFastaWriter.cpp',
  'BgzipFastqWriter.cpp',
  'BgzipWriter.cpp',
  'ChemistryTable.cpp',
  'CollectionMetadata.cpp',
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
  'ErrnoReason.cpp',
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
  'LibraryInfo.cpp',
  'MD5.cpp',
  'MemoryUtils.cpp',
  'PbiBuilder.cpp',
  'PbiBuilderBase.cpp',
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
  'SamReader.cpp',
  'SamTagCodec.cpp',
  'SamWriter.cpp',
  'SequenceInfo.cpp',
  'StringUtilities.cpp',
  'Tag.cpp',
  'TagCollection.cpp',
  'TextFileReader.cpp',
  'TextFileWriter.cpp',
  'TimeUtils.cpp',
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

  # csv
  'csv/CsvReader.cpp',
  'csv/CsvWriter.cpp',

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

pbbam_cpp_sources += pbbam_version_sources

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

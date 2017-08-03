###########
# bam2sam #
###########

pbbam_Bam2SamVersion_h_config = configuration_data()
pbbam_Bam2SamVersion_h_config.set('Bam2Sam_VERSION', meson.project_version())
pbbam_Bam2SamVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_Bam2SamVersion_h_config.set('PacBioBAM_TestsDir', join_paths(meson.source_root(), 'tests'))

configure_file(input : 'bam2sam/src/Bam2SamVersion.h.in', output : 'Bam2SamVersion.h', configuration : pbbam_Bam2SamVersion_h_config)

if pbbam_htslib_dep.found()
  pbbam_htslib_libdir = pbbam_htslib_dep.get_pkgconfig_variable('libdir')
endif

pbbam_bam2sam = executable(
  'bam2sam', [
    'common/OptionParser.cpp',
    'bam2sam/src/main.cpp',
    'bam2sam/src/Bam2Sam.cpp'
  ],
  dependencies : [pbbam_htslib_dep, pbbam_thread_dep],
  include_directories : [pbbam_include_directories, include_directories('bam2sam')],
  install : get_option('enable-pbbam-build-tools'),
  build_rpath : pbbam_htslib_libdir,
  cpp_args : pbbam_warning_flags)

# tests
if get_option('enable-pbbam-tests')
  bam2sam_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/bam2sam.t.in']),
    output : 'bam2sam.t', configuration : pbbam_Bam2SamVersion_h_config)

  test(
    'bam2sam_CramTests',
    pbbam_python,
    args : pbbam_cram_py + [bam2sam_t])
endif

###########
# pbindex #
###########

pbbam_PbIndexVersion_h_config = configuration_data()
pbbam_PbIndexVersion_h_config.set('PbIndex_VERSION', meson.project_version())
configure_file(input : 'pbindex/src/PbIndexVersion.h.in', output : 'PbIndexVersion.h', configuration : pbbam_PbIndexVersion_h_config)

pbbam_pbindex = executable(
  'pbindex', [
    'common/OptionParser.cpp',
    'pbindex/src/main.cpp',
    'pbindex/src/PbIndex.cpp'
  ],
  dependencies : pbbam_htslib_dep,
  include_directories : [pbbam_include_directories, include_directories('pbindex')],
  link_with : pbbam_lib,
  install : get_option('enable-pbbam-build-tools'),
  cpp_args : pbbam_warning_flags)

###############
# pbindexdump #
###############

pbbam_PbIndexDumpVersion_h_config = configuration_data()
pbbam_PbIndexDumpVersion_h_config.set('PbIndexDump_VERSION', meson.project_version())
pbbam_PbIndexDumpVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_PbIndexDumpVersion_h_config.set('PacBioBAM_TestsDir', join_paths(meson.source_root(), 'tests'))
configure_file(input : 'pbindexdump/src/PbIndexDumpVersion.h.in', output : 'PbIndexDumpVersion.h', configuration : pbbam_PbIndexDumpVersion_h_config)

pbbam_pbindexdump = executable(
  'pbindexdump', [
    'common/OptionParser.cpp',
    'pbindexdump/src/CppFormatter.cpp',
    'pbindexdump/src/JsonFormatter.cpp',
    'pbindexdump/src/PbIndexDump.cpp',
    'pbindexdump/src/main.cpp'
  ],
  dependencies : pbbam_htslib_dep,
  include_directories : [pbbam_include_directories, include_directories('pbindexdump')],
  link_with : pbbam_lib,
  install : get_option('enable-pbbam-build-tools'),
  cpp_args : pbbam_warning_flags)

# tests
if get_option('enable-pbbam-tests')
  pbindexdump_json_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbindexdump_json.t.in']),
    output : 'pbindexdump_json.t', configuration : pbbam_PbIndexDumpVersion_h_config)
  pbindexdump_cpp_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbindexdump_cpp.t.in']),
    output : 'pbindexdump_cpp.t', configuration : pbbam_PbIndexDumpVersion_h_config)

  test('pbindexdump_CramTests', pbbam_python, args : pbbam_cram_py + [pbindexdump_json_t, pbindexdump_cpp_t])
endif

###########
# pbmerge #
###########

pbbam_PbMergeVersion_h_config = configuration_data()
pbbam_PbMergeVersion_h_config.set('PbMerge_VERSION', meson.project_version())
pbbam_PbMergeVersion_h_config.set('PacBioBAM_VERSION', meson.project_version())
pbbam_PbMergeVersion_h_config.set('PacBioBAM_BinDir', meson.current_build_dir())
pbbam_PbMergeVersion_h_config.set('PacBioBAM_TestsDir', join_paths(meson.source_root(), 'tests'))
pbbam_PbMergeVersion_h_config.set('GeneratedTestDataDir', join_paths(meson.current_build_dir()))
configure_file(input : 'pbmerge/src/PbMergeVersion.h.in', output : 'PbMergeVersion.h', configuration : pbbam_PbMergeVersion_h_config)

pbbam_pbmerge = executable(
  'pbmerge', [
    'common/BamFileMerger.h',
    'common/OptionParser.cpp',
    'pbmerge/src/main.cpp'
  ],
  dependencies : [pbbam_boost_dep, pbbam_htslib_dep],
  include_directories : [pbbam_include_directories, include_directories('pbmerge')],
  link_with : pbbam_lib,
  install : get_option('enable-pbbam-build-tools'),
  cpp_args : pbbam_warning_flags)

# tests
if get_option('enable-pbbam-tests') and not get_option('enable-pbbam-auto-validate')
  pbmerge_pacbio_ordering_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbmerge_pacbio_ordering.t.in']),
    output : 'pbmerge_pacbio_ordering.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_aligned_ordering_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbmerge_aligned_ordering.t.in']),
    output : 'pbmerge_aligned_ordering.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_mixed_ordering_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbmerge_mixed_ordering.t.in']),
    output : 'pbmerge_mixed_ordering.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_dataset_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbmerge_dataset.t.in']),
    output : 'pbmerge_dataset.t', configuration : pbbam_PbMergeVersion_h_config)
  pbmerge_fofn_t = configure_file(
    input : join_paths([meson.source_root(), 'tests/src/cram/pbmerge_fofn.t.in']),
    output : 'pbmerge_fofn.t', configuration : pbbam_PbMergeVersion_h_config)

  test('pbmerge_CramTests', pbbam_python, args : pbbam_cram_py + [
    pbmerge_pacbio_ordering_t,
    pbmerge_aligned_ordering_t,
    pbmerge_mixed_ordering_t,
    pbmerge_dataset_t,
    pbmerge_fofn_t])
endif

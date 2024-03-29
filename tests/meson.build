subdir('src')

pbbam_python = import('python').find_installation('python3')
pbbam_cram_script = find_program('cram', required : false)
if not pbbam_cram_script.found()
  warning('Using bundled cram script')
  pbbam_cram_script = find_program('scripts/cram.py')
endif

pbbam_gtest_dep = dependency(
  'gtest',
  main : true,
  fallback : ['gtest', 'gtest_main_dep'])

pbbam_PbbamTestData_h_config = configuration_data()
pbbam_PbbamTestData_h_config.set('PacBioBAM_BinDir', join_paths([meson.current_build_dir(), '../tools']))
pbbam_PbbamTestData_h_config.set('PacBioBAM_TestsDir', meson.current_source_dir())
pbbam_PbbamTestData_h_config.set('CMAKE_CURRENT_BINARY_DIR', meson.current_build_dir())
pbbam_PbbamTestData_h_config.set('GeneratedDir', meson.current_build_dir())
pbbam_PbbamTestData_h_config.set('GeneratedTestDataDir', meson.current_build_dir())

pbbam_group_fofn_in = configure_file(
  input : files('data/group/group.fofn.in'),
  output : 'group.fofn',
  configuration : pbbam_PbbamTestData_h_config)
pbbam_PbbamTestData_h = configure_file(
  input : files('src/PbbamTestData.h.in'),
  output : 'PbbamTestData.h',
  configuration : pbbam_PbbamTestData_h_config)
pbbam_test_cpp_sources += pbbam_PbbamTestData_h

pbbam_test = executable(
  'pbbam_test',
  pbbam_test_cpp_sources,
  dependencies : [pbbam_gtest_dep, pbbam_boost_dep, pbbam_htslib_dep, pbbam_zlib_dep, pbbam_pbcopper_dep],
  include_directories : pbbam_include_directories,
  link_with : pbbam_lib,
  cpp_args : pbbam_warning_flags,
  install : false)

custom_target('pbbam_generate_data',
  output : 'input.fa',
  command : [
    pbbam_python,
    files('scripts/generate_data.py'),
    join_paths([meson.current_source_dir(), 'data']),
    meson.current_build_dir()],
  build_by_default : true,
  install : false)

pbbamify_synthetic_dataset = configure_file(
  input : files('data/pbbamify/synthetic_movie_all.subreadset.xml.in'),
  output : 'synthetic_movie_all.subreadset.xml',
  configuration : pbbam_PbbamTestData_h_config)

subdir('zero_byte_check')

#########
# tests #
#########

test(
  'pbbam formatting check',
  pbbam_clang_formatter,
  args : [
    '--all'],
  workdir : meson.project_source_root())

test(
  'pbbam gtest unittests',
  pbbam_test,
  args : [
    '--gtest_output=xml:' + join_paths(meson.project_build_root(), 'pbbam-gtest-unittests.xml')],
  env : [
    'ARGS=-V',
    'VERBOSE=1'])

test(
  'pbbam zero-byte input check',
  pbbam_cram_script,
  args : [
    '--xunit-file=' + join_paths(meson.project_build_root(), 'pbbam-cram-zero-byte-check.xml'),
    '--verbose'] + files('src/cram/zero_byte_check.t'),
  env : [
    '__ZERO_BYTE_CHECK_EXE=' + pbbam_zero_byte_check.full_path(),
  ])

pbbam_zero_byte_check_cpp_sources = files([
    'main.cpp'
])

pbbam_zero_byte_check = executable(
  'zero_byte_check',
  pbbam_zero_byte_check_cpp_sources,
  dependencies : [pbbam_boost_dep, pbbam_htslib_dep, pbbam_zlib_dep, pbbam_pbcopper_dep],
  include_directories : pbbam_include_directories,
  link_with : pbbam_lib,
  cpp_args : pbbam_warning_flags,
  install : false)

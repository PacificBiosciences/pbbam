project(
  'PacBioBAM',
  'cpp',
  version : '1.6.1',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++14',
    'b_ndebug=if-release'],
  license : 'BSD-3',
  meson_version : '>= 0.52.0')

############
# CXXFLAGS #
############

pbbam_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wduplicated-cond',
  '-Wduplicated-branches',
  '-Wlogical-op',
  '-Wrestrict',
  '-Wnull-dereference',
  '-Wuseless-cast',
  '-Wdouble-promotion',
  '-Wshadow',
  '-Wformat=1']
    if cpp.has_argument(cflag)
      pbbam_warning_flags += cflag
    endif
endforeach

pbbam_macros = []
if get_option('permissive-cigar')
  pbbam_macros += ['-DPBBAM_PERMISSIVE_CIGAR']
  warning('**********************************************')
  warning('* You have enabled allowing "M" in BAM files *')
  warning('*    This is an unsupported combination!     *')
  warning('**********************************************')
endif

################
# dependencies #
################

# threads
pbbam_thread_dep = dependency('threads', required : true)

# boost
pbbam_boost_dep = dependency('boost', include_type : 'system', required : true)

# zlib
pbbam_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])

# htslib
pbbam_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])

# pbcopper
pbbam_pbcopper_dep = dependency('pbcopper', required : true, fallback : ['pbcopper', 'pbcopper_dep'])

###########
# headers #
###########

subdir('include')

#####################
# sources + library #
#####################

subdir('src')

#########
# tests #
#########

pbbam_run_tests = (not meson.is_subproject()) and get_option('tests')
if pbbam_run_tests
  pbbam_clang_formatter = find_program('tools/check-formatting')
  pbbam_eof_check = find_program('tools/check-eof')
  subdir('tests')
endif

#########
# tools #
#########

subdir('tools')

#################
# documentation #
#################

if get_option('build-docs')
  subdir('docs')
endif

###################
# dependency info #
###################

if not meson.is_subproject()
  # need to add pbcopper into 'Requires:' field,
  # but Meson currently only allows this if it's a
  # 'pkgconfig-dependency object' and not a subproject
  pbbam_requires = []
  if pbbam_pbcopper_dep.type_name() == 'pkgconfig'
    pbbam_requires = [pbbam_pbcopper_dep]
  endif

  import('pkgconfig').generate(
    pbbam_lib,
    version : meson.project_version(),
    name : 'pbbam',
    requires : pbbam_requires,
    filebase : 'pbbam',
    description : 'Library for accessing PacBio-compatible BAM files')
endif

pbbam_dep = declare_dependency(
  include_directories : pbbam_include_directories,
  link_with : pbbam_lib,
  dependencies : [pbbam_htslib_dep, pbbam_pbcopper_dep],
  version : meson.project_version(),
  compile_args : pbbam_macros)

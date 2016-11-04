########################################################################
# CMake build script for htslib library.
########################################################################

cmake_policy(SET CMP0048 NEW)  # lets us set version in project()
project(htslib VERSION 1.1 LANGUAGES C)
cmake_minimum_required(VERSION 3.0)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
    find_package(Git QUIET REQUIRED)
    execute_process(COMMAND
        "${GIT_EXECUTABLE}" "describe" "--always" "--dirty=*"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        RESULT_VARIABLE res
        OUTPUT_VARIABLE out
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (NOT res EQUAL 0)
        message(FATAL_ERROR "Could not determine git sha1 via `git describe --always --dirty`")
    endif()
    set(LIBHTS_VERSION "${out}")
else()
    set(LIBHTS_VERSION "${htslib_VERSION}")
endif()

file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/generated/version.h
    CONTENT "#define HTS_VERSION \"${LIBHTS_VERSION}\""
)

set(LIBHTS_CFLAGS "-Wall -O2")
include(CheckCCompilerFlag)
check_c_compiler_flag("-Wno-unused-function" HAS_NO_UNUSED_FUNCTION)
if(HAS_NO_UNUSED_FUNCTION)
    set(LIBHTS_CFLAGS "${LIBHTS_CFLAGS} -Wno-unused-function")
endif()

find_package(Threads)
if (NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()

set(LIBHTS_SRCS
    htslib/kfunc.c
    htslib/knetfile.c
    htslib/kstring.c
    htslib/bgzf.c
    htslib/faidx.c
    htslib/hfile.c
    htslib/hfile_net.c
    htslib/hts.c
    htslib/regidx.c
    htslib/sam.c
    htslib/synced_bcf_reader.c
    htslib/vcf_sweep.c
    htslib/tbx.c
    htslib/vcf.c
    htslib/vcfutils.c
    htslib/cram/cram_codecs.c
    htslib/cram/cram_decode.c
    htslib/cram/cram_encode.c
    htslib/cram/cram_index.c
    htslib/cram/cram_io.c
    htslib/cram/cram_samtools.c
    htslib/cram/cram_stats.c
    htslib/cram/files.c
    htslib/cram/mFILE.c
    htslib/cram/md5.c
    htslib/cram/open_trace_file.c
    htslib/cram/pooled_alloc.c
    htslib/cram/rANS_static.c
    htslib/cram/sam_header.c
    htslib/cram/string_alloc.c
    htslib/cram/thread_pool.c
    htslib/cram/vlen.c
    htslib/cram/zfio.c
)

if (MINGW)
    list(APPEND LIBHTS_SRCS
        htslib/compat/drand48.c
        htslib/compat/fsync.c
        htslib/compat/getpagesize.c
        htslib/compat/valloc.c
    )
endif()

add_library(hts ${LIBHTS_SRCS})
#set_target_properties(hts PROPERTIES
#    ARCHIVE_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
#    RUNTIME_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
#    LIBRARY_OUTPUT_DIRECTORY ${PacBioBAM_LibDir}
#)

if (MINGW)
    set_target_properties(hts
        PROPERTIES
        COMPILE_FLAGS ${LIBHTS_CFLAGS}
        DEFINE_SYMBOL ""
        LINK_FLAGS "-Wl,--out-implib,libhts.dll.a"
        )
else()
    set_target_properties(hts
        PROPERTIES
        COMPILE_FLAGS ${LIBHTS_CFLAGS}
        DEFINE_SYMBOL ""
        )
endif()

if(WIN32)
    set(SOCKET_LIBRARIES "ws2_32")
endif()

# define symbols for projects that use htslib
set(HTSLIB_INCLUDE_DIRS
    ${CMAKE_BINARY_DIR}/generated
    ${CMAKE_CURRENT_SOURCE_DIR}/htslib
    CACHE INTERNAL
    "${PROJECT_NAME}: Include Directories"
    FORCE
)

set(HTSLIB_LIBRARIES
    $<TARGET_FILE:hts>
    CACHE INTERNAL
    "${PROJECT_NAME}: Libraries"
    FORCE
)

target_include_directories(hts
    PUBLIC
    ${HTSLIB_INCLUDE_DIRS}
    ${ZLIB_INCLUDE_DIRS}
)

target_link_libraries(hts
    PUBLIC
    ${ZLIB_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${SOCKET_LIBRARIES}
)

if(APPLE)
    # e.g. libhts.1.dylib
    set(HTSLIB_LIBRARIES_VERSIONED_LINK
        ${CMAKE_CURRENT_BINARY_DIR}/libhts.1${CMAKE_SHARED_LIBRARY_SUFFIX}
        CACHE INTERNAL
        ""
        FORCE
    )
else()
    # e.g. libhts.so.1
    set(HTSLIB_LIBRARIES_VERSIONED_LINK
        ${CMAKE_CURRENT_BINARY_DIR}/libhts${CMAKE_SHARED_LIBRARY_SUFFIX}.1
        CACHE INTERNAL
        ""
        FORCE
    )
endif()

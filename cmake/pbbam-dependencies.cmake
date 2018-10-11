
# pthreads
find_package(Threads REQUIRED)

# boost
if(NOT Boost_INCLUDE_DIRS)
    find_package(Boost REQUIRED)
endif()

# Winsock for htslib on Windows
if(WIN32)
    set(SOCKET_LIBRARIES "ws2_32")
endif()


find_package(PkgConfig REQUIRED)

# Find zlib
pkg_check_modules(ZLIB REQUIRED zlib)
message(STATUS "   zlib CFLAGS: ${ZLIB_CFLAGS}")
message(STATUS "   zlib LDFLAGS: ${ZLIB_LDFLAGS}")

# Find HTSlib
pkg_check_modules(HTSlib REQUIRED htslib)
message(STATUS "   HTSlib CFLAGS: ${HTSlib_CFLAGS}")
message(STATUS "   HTSlib LDFLAGS: ${HTSlib_LDFLAGS}")

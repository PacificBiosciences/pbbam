
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

# zlib
if(NOT ZLIB_INCLUDE_DIRS OR NOT ZLIB_LIBRARIES)
    find_package(ZLIB REQUIRED)
endif()

# htslib
if(NOT HTSLIB_INCLUDE_DIRS OR NOT HTSLIB_LIBRARIES)
    add_subdirectory(third-party/htslib external/htslib)
endif()


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
    find_package(HTSlib)
    set(hts_INCLUDE_DIRS ${HTSlib_INCLUDE_DIRS})
    set(hts_LIBRARIES    ${HTSlib_LIBRARIES})
else()    
    set(hts_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIRS})
    set(hts_LIBRARIES    ${HTSLIB_LIBRARIES})
endif()

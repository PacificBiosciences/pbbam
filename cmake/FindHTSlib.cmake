# Find HTSlib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSlib REQUIRED htslib)

# because CMake is trying to be extra clever,
# it will not properly load libraries with
# absolute paths in *_LIBRARIES
set(HTSlib_LIBRARIES "${HTSlib_LDFLAGS}")

message(STATUS "   HTSlib include dirs: ${HTSlib_INCLUDE_DIRS}")
message(STATUS "   HTSlib libraries: ${HTSlib_LIBRARIES}")

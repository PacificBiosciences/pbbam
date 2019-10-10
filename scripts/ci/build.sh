#!/usr/bin/env bash
set -vex

#########
# BUILD #
#########

# on PA, need to first build pbcopper+htslib
if [[ ${GCC_VERSION} == PA ]]; then
  pushd _deps/pbcopper
    meson \
      --default-library static \
      --libdir lib \
      --wrap-mode nofallback \
      --prefix "${bamboo_build_working_directory}/staging" \
      -Dtests=false \
      build .
    ninja -C build -v install
  popd

  wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
  tar -xjf htslib-1.9.tar.bz2
  pushd htslib-1.9
    CFLAGS="-O3" ./configure \
      --prefix="${bamboo_build_working_directory}/staging" \
      --libdir="${bamboo_build_working_directory}/staging/lib" \
      --disable-bz2 \
      --disable-gcs \
      --disable-libcurl \
      --disable-lzma \
      --disable-plugins \
      --disable-s3

    make -j install

    # clean out unneeded cruft and shared libs,
    # as -lhts will prefer shared libraries
    rm -rf ${bamboo_build_working_directory}/staging/{bin,share}
    rm -f ${bamboo_build_working_directory}/staging/lib/*.so*

    # set pkg-config variables
    export PKG_CONFIG_LIBDIR+=":${bamboo_build_working_directory}/staging/lib/pkgconfig"

    # convert `-I` to `-isystem` in pkg-config file in order not to trigger -Werror
    sed -e 's/-I/-isystem/g' -i "${bamboo_build_working_directory}/staging/lib/pkgconfig/htslib.pc"
  popd
fi

# configure
# '--wrap-mode nofallback' prevents meson from downloading
# stuff from the internet or using subprojects.
meson \
  --werror \
  --buildtype "${BUILDTYPE:-release}" \
  --default-library "${LIBRARYTYPE:-shared}" \
  --libdir lib \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --wrap-mode "${ENABLED_WRAP_MODE:-nofallback}" \
  --prefix "${PREFIX_ARG:-/}" \
  -Db_coverage="${ENABLED_COVERAGE:-false}" \
  -Db_lto="${ENABLED_LTO:-false}" \
  -Db_sanitize="${ENABLED_SANITIZERS:-none}" \
  -Dcpp_debugstl="${ENABLED_DEBUGSTL:-false}" \
  -Dtests="${ENABLED_TESTS:-false}" \
  "${CURRENT_BUILD_DIR:-build}" .

# build
ninja -C "${CURRENT_BUILD_DIR:-build}" -v

#!/usr/bin/env bash
set -vex

#########
# BUILD #
#########

# configure
# '--wrap-mode nofallback' prevents meson from downloading
# stuff from the internet or using subprojects.
meson \
  --werror \
  --buildtype "${BUILDTYPE:-release}" \
  --default-library shared \
  --libdir lib \
  --unity "${ENABLED_UNITY_BUILD:-off}" \
  --wrap-mode nofallback \
  --prefix "${PREFIX_ARG:-/usr/local}" \
  -Db_coverage="${ENABLED_COVERAGE:-false}" \
  -Db_sanitize="${ENABLED_SANITIZERS:-none}" \
  -Dtests="${ENABLED_TESTS:-false}" \
  "${CURRENT_BUILD_DIR:-build}" .

# build
ninja -C "${CURRENT_BUILD_DIR:-build}" -v

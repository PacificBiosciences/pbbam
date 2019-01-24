#!/usr/bin/env bash
set -vex

###########
# INSTALL #
###########

if [[ ${_install_image} != true ]]; then
  echo "Not installing image (branch: ${bamboo_planRepository_branchName}), returning."
  return 0
fi

if [[ ${PREFIX_ARG} ]]; then
  ## Cleaning out old installation from /mnt/software
  rm -rf "${PREFIX_ARG}"/*
fi

# *never* install with ASAN enabled
meson configure -Db_sanitize=none "${CURRENT_BUILD_DIR:-build}"

DESTDIR="${DESTDIR:-/}" ninja -C "${CURRENT_BUILD_DIR:-build}" -v install

#!/usr/bin/env bash
set -vex

################
# DEPENDENCIES #
################

## Load modules
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

module load meson
module load ninja

module load zlib
module load htslib
module load samtools

module load boost

module load cram

source scripts/ci/setup.sh

case "${bamboo_planRepository_branchName}" in
  develop|master)
    _install_image_default="${INSTALL_IMAGE:-false}"
    _create_artifact_default="${CREATE_ARTIFACT:-false}"

    export PREFIX_ARG="/mnt/software/p/pbbam/${bamboo_planRepository_branchName}"
    export BUILD_NUMBER="${bamboo_globalBuildNumber:-0}"
    ;;
esac

export _install_image="${_install_image_default:-false}"
export _create_artifact="${_create_artifact_default:-false}"


BOOST_ROOT="${BOOST_ROOT%/include}"
# unset these variables to have meson discover all
# boost-dependent variables from BOOST_ROOT alone
unset BOOST_INCLUDEDIR
unset BOOST_LIBRARYDIR

# in order to make shared libraries consumable
# by conda and other package managers
export LDFLAGS="-static-libstdc++ -static-libgcc"

source scripts/ci/build.sh
source scripts/ci/test.sh
source scripts/ci/install.sh
source scripts/ci/artifact.sh

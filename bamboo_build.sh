#!/usr/bin/env bash
set -vex

################
# DEPENDENCIES #
################

## Load modules
type module >& /dev/null || . /mnt/software/Modules/current/init/bash

module purge

[[ ${GCC_VERSION} == ICC* ]] && module use /pbi/dept/primary/modulefiles
module use /mnt/software/modulefiles

module load meson
module load ninja

module load zlib
module load htslib
module load samtools

# use the same boost version as PA (avoid ODR explosions), TAK-705
[[ ${GCC_VERSION} == ICC2017* || ${GCC_VERSION} == PA ]] && PA_BOOST_VERSION="/1.58"
module load boost${PA_BOOST_VERSION}

module load cram


export BUILD_NUMBER="0"
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
export LDFLAGS=${LDFLAGS:-"-fuse-ld=gold -static-libstdc++ -static-libgcc"}

source scripts/ci/setup.sh
source scripts/ci/build.sh
source scripts/ci/test.sh
source scripts/ci/install.sh
source scripts/ci/artifact.sh

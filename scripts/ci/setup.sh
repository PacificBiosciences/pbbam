#!/usr/bin/env bash
set -vex

export ENABLED_TESTS="true"

case "${GCC_VERSION}" in
  4.8)
    module load gtest/gcc48
    ;;

  next)
    module load gcc/8.1.0
    module load gtest
    ;;

  PA)
    # have to build htslib for PA
    module unload htslib
    module load zlib
    module load gtest/gcc48

    # load SCL GCC
    source /opt/rh/devtoolset-6/enable

    export NEXUS_PROJECT=pacbio/seq/pa/pbbam
    export NEXUS_TC=""
    ;;

  *)
    module load gcc
    module load gtest
    ;;
esac

module load ccache

export CC="ccache gcc"
export CXX="ccache g++"
export CCACHE_BASEDIR="${PWD}"

if [[ $USER == bamboo ]]; then
  export CCACHE_DIR=/mnt/secondary/Share/tmp/bamboo.${bamboo_shortPlanKey}.${bamboo_shortJobKey}.ccachedir
  export CCACHE_TEMPDIR=/scratch/bamboo.ccache_tempdir
fi

#!/usr/bin/env bash
set -vex

########
# TEST #
########

type module >& /dev/null || . /mnt/software/Modules/current/init/bash

# Note: htslib v1.7 added native long CIGAR support. pbbam "spoofs" it 
#       when running <1.7. So we'll always check the default htslib for 
#       general test success/fail, and then check pre-/post-v1.7 explicitly
#       to ensure we pass in either context (detectable at runtime).

# default htslib
ninja -C "${CURRENT_BUILD_DIR:-build}" -v test

# explicit htslib v1.6
module unload htslib
module load htslib/1.6
ninja -C "${CURRENT_BUILD_DIR:-build}" -v test

# explicit htslib v1.7
module unload htslib
module load htslib/1.7
ninja -C "${CURRENT_BUILD_DIR:-build}" -v test\

# restore default
module unload htslib
module load htslib

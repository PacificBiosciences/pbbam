  $ cat "${TESTDIR}"/../../data/zero_bytes.bam | "${__ZERO_BYTE_CHECK_EXE}" bam
  [pbbam] BAM reader ERROR: could not read from empty input: stdin

  $ cat "${TESTDIR}"/../../data/zero_bytes.sam | "${__ZERO_BYTE_CHECK_EXE}" sam
  [pbbam] SAM reader ERROR: could not read from empty input: stdin

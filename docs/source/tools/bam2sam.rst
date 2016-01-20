.. _bam2sam:

bam2sam
=======

::

  Usage: bam2sam [options] [input]

  bam2sam converts a BAM file to SAM. It is essentially a stripped-down 'samtools
  view', mostly useful for testing/debugging without requiring samtools. Input BAM
  file is read from a file or stdin, and SAM output is written to stdout.

  Options:
    -h, --help            show this help message and exit
    --version             show program's version number and exit

  Options:
    input               Input BAM file. If not provided, stdin will be used as input.
    --no-header         Omit header from output.
    --header-only       Print only the header (no records).

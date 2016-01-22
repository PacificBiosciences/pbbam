.. _pbmerge:

pbmerge
=======

::

  Usage: pbmerge [options] [-o <out.bam>] <INPUT>

  pbmerge merges PacBio BAM files. If the input is DataSetXML, any filters will be
  applied. If no output filename is specified, new BAM will be written to stdout.

  Options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

  Input/Output:
    -o output           Output BAM filename.
    --no-pbi            Set this option to skip PBI index file creation. PBI
                        creation is automatically skipped if no output filename
                        is provided.
    INPUT               Input may be one of:
                            DataSetXML, list of BAM files, or FOFN

                            fofn: pbmerge -o merged.bam bams.fofn

                            bams: pbmerge -o merged.bam 1.bam 2.bam 3.bam

                            xml:  pbmerge -o merged.bam foo.subreadset.xml


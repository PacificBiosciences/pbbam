Setup:

  $ PATH="$TESTDIR/../../../bin:$PATH" && export PATH
  $ PBMERGE="pbmerge" && export PBMERGE
  $ BAM2SAM="bam2sam" && export BAM2SAM

  $ DATADIR="$TESTDIR/../../data" && export DATADIR
  $ INPUT_XML="$DATADIR/polymerase/consolidate.subread.dataset.xml" && export INPUT_XML
  $ BAM_1="$DATADIR/polymerase/production.subreads.bam" && export BAM_1
  $ BAM_2="$DATADIR/polymerase/production.scraps.bam" && export BAM_2

  $ MERGED_BAM="/tmp/merged.bam" && export MERGED_BAM
  $ MERGED_BAM_PBI="/tmp/merged.bam.pbi" && export MERGED_BAM_PBI

Sanity Check:

  $ $BAM2SAM --no-header $BAM_1 | cut -f 1
  ArminsFakeMovie/0/2659_3025
  ArminsFakeMovie/0/3116_3628
  ArminsFakeMovie/0/3722_4267
  ArminsFakeMovie/0/4356_4864
  ArminsFakeMovie/0/4960_5477
  ArminsFakeMovie/0/5571_6087
  ArminsFakeMovie/0/6199_6719
  ArminsFakeMovie/0/6812_7034

  $ $BAM2SAM --no-header $BAM_2  | cut -f 1
  ArminsFakeMovie/0/0_2659
  ArminsFakeMovie/0/3025_3047
  ArminsFakeMovie/0/3047_3095
  ArminsFakeMovie/0/3095_3116
  ArminsFakeMovie/0/3628_3650
  ArminsFakeMovie/0/3650_3700
  ArminsFakeMovie/0/3700_3722
  ArminsFakeMovie/0/4267_4289
  ArminsFakeMovie/0/4289_4335
  ArminsFakeMovie/0/4335_4356
  ArminsFakeMovie/0/4864_4888
  ArminsFakeMovie/0/4888_4939
  ArminsFakeMovie/0/4939_4960
  ArminsFakeMovie/0/5477_5498
  ArminsFakeMovie/0/5498_5546
  ArminsFakeMovie/0/5546_5571
  ArminsFakeMovie/0/6087_6116
  ArminsFakeMovie/0/6116_6173
  ArminsFakeMovie/0/6173_6199
  ArminsFakeMovie/0/6719_6740
  ArminsFakeMovie/0/6740_6790
  ArminsFakeMovie/0/6790_6812
  ArminsFakeMovie/0/7034_7035

Normal Merge from XML:

  $ $PBMERGE -o $MERGED_BAM $INPUT_XML

  $ [ -f $MERGED_BAM ] && echo "Found" || echo "Not found"
  Found

  $ [ -f $MERGED_BAM_PBI ] && echo "Found" || echo "Not found"
  Found

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:8aaede36\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:CodecV1=ip;BINDINGKIT=FakeBindKit;SEQUENCINGKIT=FakeSeqKit;BASECALLERVERSION=0.2.0;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:CodecV1=ip;BINDINGKIT=FakeBindKit;SEQUENCINGKIT=FakeSeqKit;BASECALLERVERSION=0.2.0;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:BAZ_FORMAT\tVN:0.3.0 (esc)
  @PG\tID:PPA-BAZ2BAM\tVN:0.1.0 (esc)
  @PG\tID:PPA-BAZWRITER\tVN:0.2.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/0/4267_4289
  ArminsFakeMovie/0/4289_4335
  ArminsFakeMovie/0/4335_4356
  ArminsFakeMovie/0/4356_4864
  ArminsFakeMovie/0/4864_4888
  ArminsFakeMovie/0/4888_4939
  ArminsFakeMovie/0/4939_4960
  ArminsFakeMovie/0/4960_5477

  $ rm $MERGED_BAM
  $ rm $MERGED_BAM_PBI

Normal Merge from XML (disabled PBI):

  $ $PBMERGE --no-pbi -o $MERGED_BAM $INPUT_XML

  $ [ -f $MERGED_BAM ] && echo "Found" || echo "Not found"
  Found

  $ [ -f $MERGED_BAM_PBI ] && echo "Found" || echo "Not found"
  Not found

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:8aaede36\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:CodecV1=ip;BINDINGKIT=FakeBindKit;SEQUENCINGKIT=FakeSeqKit;BASECALLERVERSION=0.2.0;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:CodecV1=ip;BINDINGKIT=FakeBindKit;SEQUENCINGKIT=FakeSeqKit;BASECALLERVERSION=0.2.0;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:BAZ_FORMAT\tVN:0.3.0 (esc)
  @PG\tID:PPA-BAZ2BAM\tVN:0.1.0 (esc)
  @PG\tID:PPA-BAZWRITER\tVN:0.2.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/0/4267_4289
  ArminsFakeMovie/0/4289_4335
  ArminsFakeMovie/0/4335_4356
  ArminsFakeMovie/0/4356_4864
  ArminsFakeMovie/0/4864_4888
  ArminsFakeMovie/0/4888_4939
  ArminsFakeMovie/0/4939_4960
  ArminsFakeMovie/0/4960_5477

  $ rm $MERGED_BAM

Write to stdout:

  $ $PBMERGE --no-pbi $INPUT_XML > $MERGED_BAM

  $ [ -f $MERGED_BAM ] && echo "Found" || echo "Not found"
  Found

  $ [ -f $MERGED_BAM_PBI ] && echo "Found" || echo "Not found"
  Not found

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:8aaede36\tPL:PACBIO\tDS:READTYPE=SUBREAD;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:CodecV1=ip;BINDINGKIT=FakeBindKit;SEQUENCINGKIT=FakeSeqKit;BASECALLERVERSION=0.2.0;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:CodecV1=ip;BINDINGKIT=FakeBindKit;SEQUENCINGKIT=FakeSeqKit;BASECALLERVERSION=0.2.0;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:BAZ_FORMAT\tVN:0.3.0 (esc)
  @PG\tID:PPA-BAZ2BAM\tVN:0.1.0 (esc)
  @PG\tID:PPA-BAZWRITER\tVN:0.2.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/0/4267_4289
  ArminsFakeMovie/0/4289_4335
  ArminsFakeMovie/0/4335_4356
  ArminsFakeMovie/0/4356_4864
  ArminsFakeMovie/0/4864_4888
  ArminsFakeMovie/0/4888_4939
  ArminsFakeMovie/0/4939_4960
  ArminsFakeMovie/0/4960_5477

  $ rm $MERGED_BAM

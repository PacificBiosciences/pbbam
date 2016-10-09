Setup:

  $ PATH="$TESTDIR/../../../bin:$PATH" && export PATH
  $ PBMERGE="pbmerge" && export PBMERGE
  $ BAM2SAM="bam2sam" && export BAM2SAM

  $ DATADIR="$TESTDIR/../../data" && export DATADIR
  $ HQREGION_BAM="$DATADIR/polymerase/internal.hqregions.bam" && export HQREGION_BAM
  $ SCRAPS_BAM="$DATADIR/polymerase/internal.scraps.bam" && export SCRAPS_BAM

  $ MERGED_BAM="/tmp/pacbio_ordering_merged.bam" && export MERGED_BAM
  $ MERGED_BAM_PBI="/tmp/pacbio_ordering_merged.bam.pbi" && export MERGED_BAM_PBI

Sanity Check:

  $ $BAM2SAM --header-only $HQREGION_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:ca75d884\tPL:PACBIO\tDS:READTYPE=HQREGION;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)

  $ $BAM2SAM --no-header $HQREGION_BAM | cut -f 1
  ArminsFakeMovie/100000/2659_7034

  $ $BAM2SAM --header-only $SCRAPS_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)

  $ $BAM2SAM --no-header $SCRAPS_BAM | cut -f 1
  ArminsFakeMovie/100000/0_2659
  ArminsFakeMovie/100000/3025_3047
  ArminsFakeMovie/100000/3047_3095
  ArminsFakeMovie/100000/3095_3116
  ArminsFakeMovie/100000/3628_3650
  ArminsFakeMovie/100000/3650_3700
  ArminsFakeMovie/100000/3700_3722
  ArminsFakeMovie/100000/4267_4289
  ArminsFakeMovie/100000/4289_4335
  ArminsFakeMovie/100000/4335_4356
  ArminsFakeMovie/100000/4864_4888
  ArminsFakeMovie/100000/4888_4939
  ArminsFakeMovie/100000/4939_4960
  ArminsFakeMovie/100000/5477_5498
  ArminsFakeMovie/100000/5498_5546
  ArminsFakeMovie/100000/5546_5571
  ArminsFakeMovie/100000/6087_6116
  ArminsFakeMovie/100000/6116_6173
  ArminsFakeMovie/100000/6173_6199
  ArminsFakeMovie/100000/6719_6740
  ArminsFakeMovie/100000/6740_6790
  ArminsFakeMovie/100000/6790_6812
  ArminsFakeMovie/100000/7034_7035
  ArminsFakeMovie/200000/0_2659
  ArminsFakeMovie/200000/3025_3047
  ArminsFakeMovie/200000/3047_3095
  ArminsFakeMovie/200000/3095_3116
  ArminsFakeMovie/200000/3628_3650
  ArminsFakeMovie/200000/3650_3700
  ArminsFakeMovie/200000/3700_3722
  ArminsFakeMovie/200000/4267_4289
  ArminsFakeMovie/200000/4289_4335
  ArminsFakeMovie/200000/4335_4356
  ArminsFakeMovie/200000/4864_4888
  ArminsFakeMovie/200000/4888_4939
  ArminsFakeMovie/200000/4939_4960
  ArminsFakeMovie/200000/5477_5498
  ArminsFakeMovie/200000/5498_5546
  ArminsFakeMovie/200000/5546_5571
  ArminsFakeMovie/200000/6087_6116
  ArminsFakeMovie/200000/6116_6173
  ArminsFakeMovie/200000/6173_6199
  ArminsFakeMovie/200000/6719_6740
  ArminsFakeMovie/200000/6740_6790
  ArminsFakeMovie/200000/6790_6812
  ArminsFakeMovie/200000/7034_7035
  ArminsFakeMovie/300000/0_2659
  ArminsFakeMovie/300000/3025_3047
  ArminsFakeMovie/300000/3047_3095
  ArminsFakeMovie/300000/3095_3116
  ArminsFakeMovie/300000/3628_3650
  ArminsFakeMovie/300000/3650_3700
  ArminsFakeMovie/300000/3700_3722
  ArminsFakeMovie/300000/4267_4289
  ArminsFakeMovie/300000/4289_4335
  ArminsFakeMovie/300000/4335_4356
  ArminsFakeMovie/300000/4864_4888
  ArminsFakeMovie/300000/4888_4939
  ArminsFakeMovie/300000/4939_4960
  ArminsFakeMovie/300000/5477_5498
  ArminsFakeMovie/300000/5498_5546
  ArminsFakeMovie/300000/5546_5571
  ArminsFakeMovie/300000/6087_6116
  ArminsFakeMovie/300000/6116_6173
  ArminsFakeMovie/300000/6173_6199
  ArminsFakeMovie/300000/6719_6740
  ArminsFakeMovie/300000/6740_6790
  ArminsFakeMovie/300000/6790_6812
  ArminsFakeMovie/300000/7034_7035

Normal Merge:

  $ $PBMERGE $HQREGION_BAM $SCRAPS_BAM > $MERGED_BAM

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:ca75d884\tPL:PACBIO\tDS:READTYPE=HQREGION;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/100000/0_2659
  ArminsFakeMovie/100000/2659_7034
  ArminsFakeMovie/100000/3025_3047
  ArminsFakeMovie/100000/3047_3095
  ArminsFakeMovie/100000/3095_3116
  ArminsFakeMovie/100000/3628_3650
  ArminsFakeMovie/100000/3650_3700
  ArminsFakeMovie/100000/3700_3722
  ArminsFakeMovie/100000/4267_4289
  ArminsFakeMovie/100000/4289_4335
  ArminsFakeMovie/100000/4335_4356
  ArminsFakeMovie/100000/4864_4888
  ArminsFakeMovie/100000/4888_4939
  ArminsFakeMovie/100000/4939_4960
  ArminsFakeMovie/100000/5477_5498
  ArminsFakeMovie/100000/5498_5546
  ArminsFakeMovie/100000/5546_5571
  ArminsFakeMovie/100000/6087_6116
  ArminsFakeMovie/100000/6116_6173
  ArminsFakeMovie/100000/6173_6199
  ArminsFakeMovie/100000/6719_6740
  ArminsFakeMovie/100000/6740_6790
  ArminsFakeMovie/100000/6790_6812
  ArminsFakeMovie/100000/7034_7035
  ArminsFakeMovie/200000/0_2659
  ArminsFakeMovie/200000/3025_3047
  ArminsFakeMovie/200000/3047_3095
  ArminsFakeMovie/200000/3095_3116
  ArminsFakeMovie/200000/3628_3650
  ArminsFakeMovie/200000/3650_3700
  ArminsFakeMovie/200000/3700_3722
  ArminsFakeMovie/200000/4267_4289
  ArminsFakeMovie/200000/4289_4335
  ArminsFakeMovie/200000/4335_4356
  ArminsFakeMovie/200000/4864_4888
  ArminsFakeMovie/200000/4888_4939
  ArminsFakeMovie/200000/4939_4960
  ArminsFakeMovie/200000/5477_5498
  ArminsFakeMovie/200000/5498_5546
  ArminsFakeMovie/200000/5546_5571
  ArminsFakeMovie/200000/6087_6116
  ArminsFakeMovie/200000/6116_6173
  ArminsFakeMovie/200000/6173_6199
  ArminsFakeMovie/200000/6719_6740
  ArminsFakeMovie/200000/6740_6790
  ArminsFakeMovie/200000/6790_6812
  ArminsFakeMovie/200000/7034_7035
  ArminsFakeMovie/300000/0_2659
  ArminsFakeMovie/300000/3025_3047
  ArminsFakeMovie/300000/3047_3095
  ArminsFakeMovie/300000/3095_3116
  ArminsFakeMovie/300000/3628_3650
  ArminsFakeMovie/300000/3650_3700
  ArminsFakeMovie/300000/3700_3722
  ArminsFakeMovie/300000/4267_4289
  ArminsFakeMovie/300000/4289_4335
  ArminsFakeMovie/300000/4335_4356
  ArminsFakeMovie/300000/4864_4888
  ArminsFakeMovie/300000/4888_4939
  ArminsFakeMovie/300000/4939_4960
  ArminsFakeMovie/300000/5477_5498
  ArminsFakeMovie/300000/5498_5546
  ArminsFakeMovie/300000/5546_5571
  ArminsFakeMovie/300000/6087_6116
  ArminsFakeMovie/300000/6116_6173
  ArminsFakeMovie/300000/6173_6199
  ArminsFakeMovie/300000/6719_6740
  ArminsFakeMovie/300000/6740_6790
  ArminsFakeMovie/300000/6790_6812
  ArminsFakeMovie/300000/7034_7035

  $ rm $MERGED_BAM

Shuffle Input:

  $ $PBMERGE $SCRAPS_BAM $HQREGION_BAM  > $MERGED_BAM

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:ca75d884\tPL:PACBIO\tDS:READTYPE=HQREGION;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/100000/0_2659
  ArminsFakeMovie/100000/2659_7034
  ArminsFakeMovie/100000/3025_3047
  ArminsFakeMovie/100000/3047_3095
  ArminsFakeMovie/100000/3095_3116
  ArminsFakeMovie/100000/3628_3650
  ArminsFakeMovie/100000/3650_3700
  ArminsFakeMovie/100000/3700_3722
  ArminsFakeMovie/100000/4267_4289
  ArminsFakeMovie/100000/4289_4335
  ArminsFakeMovie/100000/4335_4356
  ArminsFakeMovie/100000/4864_4888
  ArminsFakeMovie/100000/4888_4939
  ArminsFakeMovie/100000/4939_4960
  ArminsFakeMovie/100000/5477_5498
  ArminsFakeMovie/100000/5498_5546
  ArminsFakeMovie/100000/5546_5571
  ArminsFakeMovie/100000/6087_6116
  ArminsFakeMovie/100000/6116_6173
  ArminsFakeMovie/100000/6173_6199
  ArminsFakeMovie/100000/6719_6740
  ArminsFakeMovie/100000/6740_6790
  ArminsFakeMovie/100000/6790_6812
  ArminsFakeMovie/100000/7034_7035
  ArminsFakeMovie/200000/0_2659
  ArminsFakeMovie/200000/3025_3047
  ArminsFakeMovie/200000/3047_3095
  ArminsFakeMovie/200000/3095_3116
  ArminsFakeMovie/200000/3628_3650
  ArminsFakeMovie/200000/3650_3700
  ArminsFakeMovie/200000/3700_3722
  ArminsFakeMovie/200000/4267_4289
  ArminsFakeMovie/200000/4289_4335
  ArminsFakeMovie/200000/4335_4356
  ArminsFakeMovie/200000/4864_4888
  ArminsFakeMovie/200000/4888_4939
  ArminsFakeMovie/200000/4939_4960
  ArminsFakeMovie/200000/5477_5498
  ArminsFakeMovie/200000/5498_5546
  ArminsFakeMovie/200000/5546_5571
  ArminsFakeMovie/200000/6087_6116
  ArminsFakeMovie/200000/6116_6173
  ArminsFakeMovie/200000/6173_6199
  ArminsFakeMovie/200000/6719_6740
  ArminsFakeMovie/200000/6740_6790
  ArminsFakeMovie/200000/6790_6812
  ArminsFakeMovie/200000/7034_7035
  ArminsFakeMovie/300000/0_2659
  ArminsFakeMovie/300000/3025_3047
  ArminsFakeMovie/300000/3047_3095
  ArminsFakeMovie/300000/3095_3116
  ArminsFakeMovie/300000/3628_3650
  ArminsFakeMovie/300000/3650_3700
  ArminsFakeMovie/300000/3700_3722
  ArminsFakeMovie/300000/4267_4289
  ArminsFakeMovie/300000/4289_4335
  ArminsFakeMovie/300000/4335_4356
  ArminsFakeMovie/300000/4864_4888
  ArminsFakeMovie/300000/4888_4939
  ArminsFakeMovie/300000/4939_4960
  ArminsFakeMovie/300000/5477_5498
  ArminsFakeMovie/300000/5498_5546
  ArminsFakeMovie/300000/5546_5571
  ArminsFakeMovie/300000/6087_6116
  ArminsFakeMovie/300000/6116_6173
  ArminsFakeMovie/300000/6173_6199
  ArminsFakeMovie/300000/6719_6740
  ArminsFakeMovie/300000/6740_6790
  ArminsFakeMovie/300000/6790_6812
  ArminsFakeMovie/300000/7034_7035

  $ rm $MERGED_BAM

Explicit Output Filename (also enables PBI):

  $ $PBMERGE -o $MERGED_BAM $HQREGION_BAM $SCRAPS_BAM

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:ca75d884\tPL:PACBIO\tDS:READTYPE=HQREGION;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/100000/0_2659
  ArminsFakeMovie/100000/2659_7034
  ArminsFakeMovie/100000/3025_3047
  ArminsFakeMovie/100000/3047_3095
  ArminsFakeMovie/100000/3095_3116
  ArminsFakeMovie/100000/3628_3650
  ArminsFakeMovie/100000/3650_3700
  ArminsFakeMovie/100000/3700_3722
  ArminsFakeMovie/100000/4267_4289
  ArminsFakeMovie/100000/4289_4335
  ArminsFakeMovie/100000/4335_4356
  ArminsFakeMovie/100000/4864_4888
  ArminsFakeMovie/100000/4888_4939
  ArminsFakeMovie/100000/4939_4960
  ArminsFakeMovie/100000/5477_5498
  ArminsFakeMovie/100000/5498_5546
  ArminsFakeMovie/100000/5546_5571
  ArminsFakeMovie/100000/6087_6116
  ArminsFakeMovie/100000/6116_6173
  ArminsFakeMovie/100000/6173_6199
  ArminsFakeMovie/100000/6719_6740
  ArminsFakeMovie/100000/6740_6790
  ArminsFakeMovie/100000/6790_6812
  ArminsFakeMovie/100000/7034_7035
  ArminsFakeMovie/200000/0_2659
  ArminsFakeMovie/200000/3025_3047
  ArminsFakeMovie/200000/3047_3095
  ArminsFakeMovie/200000/3095_3116
  ArminsFakeMovie/200000/3628_3650
  ArminsFakeMovie/200000/3650_3700
  ArminsFakeMovie/200000/3700_3722
  ArminsFakeMovie/200000/4267_4289
  ArminsFakeMovie/200000/4289_4335
  ArminsFakeMovie/200000/4335_4356
  ArminsFakeMovie/200000/4864_4888
  ArminsFakeMovie/200000/4888_4939
  ArminsFakeMovie/200000/4939_4960
  ArminsFakeMovie/200000/5477_5498
  ArminsFakeMovie/200000/5498_5546
  ArminsFakeMovie/200000/5546_5571
  ArminsFakeMovie/200000/6087_6116
  ArminsFakeMovie/200000/6116_6173
  ArminsFakeMovie/200000/6173_6199
  ArminsFakeMovie/200000/6719_6740
  ArminsFakeMovie/200000/6740_6790
  ArminsFakeMovie/200000/6790_6812
  ArminsFakeMovie/200000/7034_7035
  ArminsFakeMovie/300000/0_2659
  ArminsFakeMovie/300000/3025_3047
  ArminsFakeMovie/300000/3047_3095
  ArminsFakeMovie/300000/3095_3116
  ArminsFakeMovie/300000/3628_3650
  ArminsFakeMovie/300000/3650_3700
  ArminsFakeMovie/300000/3700_3722
  ArminsFakeMovie/300000/4267_4289
  ArminsFakeMovie/300000/4289_4335
  ArminsFakeMovie/300000/4335_4356
  ArminsFakeMovie/300000/4864_4888
  ArminsFakeMovie/300000/4888_4939
  ArminsFakeMovie/300000/4939_4960
  ArminsFakeMovie/300000/5477_5498
  ArminsFakeMovie/300000/5498_5546
  ArminsFakeMovie/300000/5546_5571
  ArminsFakeMovie/300000/6087_6116
  ArminsFakeMovie/300000/6116_6173
  ArminsFakeMovie/300000/6173_6199
  ArminsFakeMovie/300000/6719_6740
  ArminsFakeMovie/300000/6740_6790
  ArminsFakeMovie/300000/6790_6812
  ArminsFakeMovie/300000/7034_7035

  $ [ -f $MERGED_BAM_PBI ] && echo "Found" || echo "Not found"
  Found

  $ rm $MERGED_BAM
  $ rm $MERGED_BAM_PBI

Explicit Output Filename (with disabled PBI):

  $ $PBMERGE -o $MERGED_BAM --no-pbi $HQREGION_BAM $SCRAPS_BAM

  $ $BAM2SAM --header-only $MERGED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:ca75d884\tPL:PACBIO\tDS:READTYPE=HQREGION;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @RG\tID:e83fc9c6\tPL:PACBIO\tDS:READTYPE=SCRAP;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie\tPM:SEQUEL (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)
  @PG\tID:pbmerge-0.7.0\tPN:pbmerge\tVN:0.7.0 (esc)

  $ $BAM2SAM --no-header $MERGED_BAM | cut -f 1
  ArminsFakeMovie/100000/0_2659
  ArminsFakeMovie/100000/2659_7034
  ArminsFakeMovie/100000/3025_3047
  ArminsFakeMovie/100000/3047_3095
  ArminsFakeMovie/100000/3095_3116
  ArminsFakeMovie/100000/3628_3650
  ArminsFakeMovie/100000/3650_3700
  ArminsFakeMovie/100000/3700_3722
  ArminsFakeMovie/100000/4267_4289
  ArminsFakeMovie/100000/4289_4335
  ArminsFakeMovie/100000/4335_4356
  ArminsFakeMovie/100000/4864_4888
  ArminsFakeMovie/100000/4888_4939
  ArminsFakeMovie/100000/4939_4960
  ArminsFakeMovie/100000/5477_5498
  ArminsFakeMovie/100000/5498_5546
  ArminsFakeMovie/100000/5546_5571
  ArminsFakeMovie/100000/6087_6116
  ArminsFakeMovie/100000/6116_6173
  ArminsFakeMovie/100000/6173_6199
  ArminsFakeMovie/100000/6719_6740
  ArminsFakeMovie/100000/6740_6790
  ArminsFakeMovie/100000/6790_6812
  ArminsFakeMovie/100000/7034_7035
  ArminsFakeMovie/200000/0_2659
  ArminsFakeMovie/200000/3025_3047
  ArminsFakeMovie/200000/3047_3095
  ArminsFakeMovie/200000/3095_3116
  ArminsFakeMovie/200000/3628_3650
  ArminsFakeMovie/200000/3650_3700
  ArminsFakeMovie/200000/3700_3722
  ArminsFakeMovie/200000/4267_4289
  ArminsFakeMovie/200000/4289_4335
  ArminsFakeMovie/200000/4335_4356
  ArminsFakeMovie/200000/4864_4888
  ArminsFakeMovie/200000/4888_4939
  ArminsFakeMovie/200000/4939_4960
  ArminsFakeMovie/200000/5477_5498
  ArminsFakeMovie/200000/5498_5546
  ArminsFakeMovie/200000/5546_5571
  ArminsFakeMovie/200000/6087_6116
  ArminsFakeMovie/200000/6116_6173
  ArminsFakeMovie/200000/6173_6199
  ArminsFakeMovie/200000/6719_6740
  ArminsFakeMovie/200000/6740_6790
  ArminsFakeMovie/200000/6790_6812
  ArminsFakeMovie/200000/7034_7035
  ArminsFakeMovie/300000/0_2659
  ArminsFakeMovie/300000/3025_3047
  ArminsFakeMovie/300000/3047_3095
  ArminsFakeMovie/300000/3095_3116
  ArminsFakeMovie/300000/3628_3650
  ArminsFakeMovie/300000/3650_3700
  ArminsFakeMovie/300000/3700_3722
  ArminsFakeMovie/300000/4267_4289
  ArminsFakeMovie/300000/4289_4335
  ArminsFakeMovie/300000/4335_4356
  ArminsFakeMovie/300000/4864_4888
  ArminsFakeMovie/300000/4888_4939
  ArminsFakeMovie/300000/4939_4960
  ArminsFakeMovie/300000/5477_5498
  ArminsFakeMovie/300000/5498_5546
  ArminsFakeMovie/300000/5546_5571
  ArminsFakeMovie/300000/6087_6116
  ArminsFakeMovie/300000/6116_6173
  ArminsFakeMovie/300000/6173_6199
  ArminsFakeMovie/300000/6719_6740
  ArminsFakeMovie/300000/6740_6790
  ArminsFakeMovie/300000/6790_6812
  ArminsFakeMovie/300000/7034_7035

  $ [ -f $MERGED_BAM_PBI ] && echo "Found" || echo "Not found"
  Not found

  $ rm $MERGED_BAM

Setup:

  $ PATH="$TESTDIR/../../../bin:$PATH" && export PATH
  $ PBMERGE="pbmerge" && export PBMERGE
  $ BAM2SAM="bam2sam" && export BAM2SAM

  $ DATADIR="$TESTDIR/../../data" && export DATADIR
  $ UNALIGNED_BAM="$DATADIR/polymerase/internal.hqregions.bam" && export UNALIGNED_BAM
  $ ALIGNED_BAM="$DATADIR/dataset/bam_mapping_1.bam" && export ALIGNED_BAM

  $ MERGED_BAM="/tmp/mixed_ordering_merged.bam" && export MERGED_BAM

Sanity Check:

  $ $BAM2SAM --header-only $UNALIGNED_BAM
  @HD\tVN:1.1\tSO:unknown\tpb:3.0.1 (esc)
  @RG\tID:ca75d884\tPL:PACBIO\tDS:READTYPE=HQREGION;DeletionQV=dq;DeletionTag=dt;InsertionQV=iq;MergeQV=mq;SubstitutionQV=sq;SubstitutionTag=st;Ipd:Frames=ip;PulseWidth:Frames=pw;PkMid=pm;PkMean=pa;LabelQV=pq;AltLabel=pt;AltLabelQV=pv;PulseMergeQV=pg;PulseCall=pc;PrePulseFrames=pd;PulseCallWidth=px;BINDINGKIT=100372700;SEQUENCINGKIT=100356200;BASECALLERVERSION=0.1;FRAMERATEHZ=100.000000\tPU:ArminsFakeMovie (esc)
  @PG\tID:baz2bam-0.15.0\tPN:baz2bam\tVN:0.15.0 (esc)
  @PG\tID:bazFormat-0.3.0\tPN:bazFormat\tVN:0.3.0 (esc)
  @PG\tID:bazwriter-0.15.0\tPN:bazwriter\tVN:0.15.0 (esc)

  $ $BAM2SAM --no-header $UNALIGNED_BAM | cut -f 1
  ArminsFakeMovie/100000/2659_7034

  $ $BAM2SAM --header-only $ALIGNED_BAM
  @HD\tVN:1.3.1\tSO:coordinate\tpb:3.0.1 (esc)
  @SQ\tSN:lambda_NEB3011\tLN:48502\tM5:a1319ff90e994c8190a4fe6569d0822a (esc)
  @RG\tID:a9a22406c5\tDS:READTYPE=SUBREAD;BINDINGKIT=100356300;SEQUENCINGKIT=100356200;BASECALLERVERSION=2.3;InsertionQV=iq;DeletionQV=dq;SubstitutionQV=sq;MergeQV=mq;SubstitutionTag=st;DeletionTag=dt\tPL:PACBIO\tPU:m140905_042212_sidney_c100564852550000001823085912221377_s1_X0\tSM:c100564852550000001823085912221377 (esc)
  @PG\tID:BLASR\tVN:1.3.1.141565\tCL:/home/UNIXHOME/yli/for_the_people/blasr_bam_out/blasr m140905_042212_sidney_c100564852550000001823085912221377_s1_X0.1.bax.h5 lambdaNEB.fa -out tmp.bam -bam -bestn 10 -minMatch 12 -nproc 8 -minSubreadLength 50 -minReadLength 50 -randomSeed 1 -clipping subread  (esc)

  $ $BAM2SAM --no-header $ALIGNED_BAM | cut -f 1,3,4 | head -n 10
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/49050/48_1132\tlambda_NEB3011\t1 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/32328/0_344\tlambda_NEB3011\t676 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/6469/9936_10187\tlambda_NEB3011\t2171 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/6469/10232_10394\tlambda_NEB3011\t2204 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/30983/7468_8906\tlambda_NEB3011\t3573 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/13473/5557_7235\tlambda_NEB3011\t4507 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/13473/7285_8657\tlambda_NEB3011\t4508 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/19915/426_1045\tlambda_NEB3011\t4593 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/30983/7064_7421\tlambda_NEB3011\t4670 (esc)
  m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/19915/0_382\tlambda_NEB3011\t4843 (esc)

Normal Merge - should fail:

  $ $PBMERGE $UNALIGNED_BAM $ALIGNED_BAM > $MERGED_BAM
  ERROR: BAM file sort orders do not match, aborting merge
  [1]

Shuffle Input - should fail:

  $ $PBMERGE $ALIGNED_BAM $UNALIGNED_BAM > $MERGED_BAM
  ERROR: BAM file sort orders do not match, aborting merge
  [1]

Cleanup:

  $ rm $MERGED_BAM

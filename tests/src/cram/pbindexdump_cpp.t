Setup:

  $ PATH="$TESTDIR/../../../bin:$PATH" && export PATH

  $ PBINDEXDUMP="pbindexdump" && export PBINDEXDUMP

  $ DATADIR="$TESTDIR/../../data" && export DATADIR

Normal C++:

  $ $PBINDEXDUMP --format=cpp $DATADIR/polymerase/production_hq.hqregion.bam.pbi
  PbiRawData rawData;
  rawData.Version(PbiFile::Version_3_0_1);
  rawData.FileSections(PbiFile::BASIC);
  rawData.NumReads(1);
  
  PbiRawBasicData& basicData = rawData.BasicData();
  basicData.rgId_       = {-898246524};
  basicData.qStart_     = {2659};
  basicData.qEnd_       = {7034};
  basicData.holeNumber_ = {0};
  basicData.readQual_   = {0.01};
  basicData.ctxtFlag_   = {0};
  basicData.fileOffset_ = {20054016};
  
  
--(leave the blank lines above this)--

Request C++, with JSON options (stdout includes usage/help, so we just want to check stderr):

  $ $PBINDEXDUMP --format=cpp --json-indent-level=2 $DATADIR/polymerase/production_hq.hqregion.bam.pbi > /dev/null
  
  ERROR: JSON formatting options not valid on non-JSON output
  
  [1]

  $ $PBINDEXDUMP --format=cpp --json-raw $DATADIR/polymerase/production_hq.hqregion.bam.pbi > /dev/null
  
  ERROR: JSON formatting options not valid on non-JSON output
  
  [1]

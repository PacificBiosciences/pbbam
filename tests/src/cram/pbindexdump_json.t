Setup:

  $ PATH="$TESTDIR/../../../bin:$PATH" && export PATH

  $ PBINDEXDUMP="pbindexdump" && export PBINDEXDUMP

  $ DATADIR="$TESTDIR/../../data" && export DATADIR

Default settings (JSON):

  $ $PBINDEXDUMP $DATADIR/polymerase/production_hq.hqregion.bam.pbi
  {
      "fileSections": [
          "BasicData"
      ],
      "numReads": 1,
      "reads": [
          {
              "contextFlag": 0,
              "fileOffset": 20054016,
              "holeNumber": 0,
              "qEnd": 7034,
              "qStart": 2659,
              "readQuality": 0.00999999977648258,
              "rgId": -898246524
          }
      ],
      "version": "3.0.1"
  }

JSON indent level(2):

  $ $PBINDEXDUMP --json-indent-level=2 $DATADIR/polymerase/production_hq.hqregion.bam.pbi
  {
    "fileSections": [
      "BasicData"
    ],
    "numReads": 1,
    "reads": [
      {
        "contextFlag": 0,
        "fileOffset": 20054016,
        "holeNumber": 0,
        "qEnd": 7034,
        "qStart": 2659,
        "readQuality": 0.00999999977648258,
        "rgId": -898246524
      }
    ],
    "version": "3.0.1"
  }

JSON raw:

  $ $PBINDEXDUMP --json-raw $DATADIR/polymerase/production_hq.hqregion.bam.pbi
  {
      "basicData": {
          "ctxtFlag": [
              0
          ],
          "fileOffset": [
              20054016
          ],
          "holeNumber": [
              0
          ],
          "qEnd": [
              7034
          ],
          "qStart": [
              2659
          ],
          "readQual": [
              0.00999999977648258
          ],
          "rgId": [
              -898246524
          ]
      },
      "fileSections": [
          "BasicData"
      ],
      "numReads": 1,
      "version": "3.0.1"
  }

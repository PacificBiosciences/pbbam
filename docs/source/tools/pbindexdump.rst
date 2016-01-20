.. _pbindexdump:

pbindexdump
===========

::

  Usage: pbindexdump [options] [input]

  pbindexdump prints a human-readable view of PBI data to stdout.

  Options:
    -h, --help            show this help message and exit
    --version             show program's version number and exit

  Input/Output:
    input               Input PBI file. If not provided, stdin will be used as input.
    --format=STRING     Output format, one of:
                            json, cpp

                        json: pretty-printed JSON [default]

                        cpp: copy/paste-able C++ code that can be used to
                        construct the equivalent PacBio::BAM::PbiRawData object

  JSON Formatting:
    --json-indent-level=INT
                        JSON indent level [4]
    --json-raw          Prints fields in a manner that more closely reflects the
                        PBI file format - presenting data as per-field columns,
                        not per-record objects.

JSON Output Schemas
-------------------

Normal JSON:

.. code-block:: JSON

    {
      "type": "object",
      "properties": {
        "fileSections": {
          "type": "array",
          "items": { "type": "string" },
        },
        "numReads": { "type": "integer" },
        "reads": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "aEnd": { "type": "integer" },
              "aStart": { "type": "integer" },
              "bcForward": { "type": "integer" },
              "bcQuality": { "type": "integer" },
              "bcReverse": { "type": "integer" },
              "contextFlag": { "type": "integer" },
              "fileOffset": { "type": "integer" },
              "holeNumber": { "type": "integer" },
              "mapQuality": { "type": "integer" },
              "nM": { "type": "integer" },
              "nMM": { "type": "integer" },
              "qEnd": { "type": "integer" },
              "qStart": { "type": "integer" },
              "readQuality": { "type": "number" },
              "reverseStrand": { "type": "integer" },
              "rgId": { "type": "integer" },
              "tEnd": { "type": "integer" },
              "tId": { "type": "integer" },
              "tStart: { "type": "integer" }
            },
            "required": [
              "contextFlag",
              "fileOffset",
              "holeNumber",
              "qEnd",
              "qStart",
              "readQuality",
              "rgId"
            ]
          }
        },
        "references": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "beginRow": { "type": "integer" },
              "endRow": { "type": "integer" },
              "tId": { "type": "integer" }
            },
            "required" : [ "beginRow", "endRow","tId" ]
          }
        }q
        "version": { "type": "string" }
      },
      "required": [
        "fileSections",
        "numReads",
        "reads",
        "version"
      ]
    }

"Raw" JSON:

.. code-block:: JSON

    {
      "type": "object",
      "properties": {
        "barcodeData" : {
          "type" : "object",
          "properties: {
            "bcForward" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "bcQuality" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "bcReverse" : {
              "type": "array",
              "items" : { "type": "integer" }
            }
          }
        },
        "basicData" : {
          "type" : "object",
          "properties: {
            "contextFlag" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "fileOffset" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "holeNumber" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "qEnd" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "qStart" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "readQuality" : {
              "type": "array",
              "items" : { "type": "number" }
            },
            "rgId : {
              "type": "array",
              "items" : { "type": "integer" }
            }
          }
        },
        "fileSections": {
          "type": "array",
          "items": { "type": "string" },
        },
        "mappedData" : {
          "type" : "object",
          "properties: {
            "aEnd" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "aStart" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "mapQuality" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "nM" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "nMM" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "readQuality" : {
              "type": "array",
              "items" : { "type": "number" }
            },
            "reverseStrand" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "tEnd" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "tId" : {
              "type": "array",
              "items" : { "type": "integer" }
            },
            "tStart" : {
              "type": "array",
              "items" : { "type": "integer" }
            }
          }
        },
        "numReads": { "type": "integer" },
        "references": {
          "type": "array",
          "items": {
            "type": "object",
            "properties": {
              "beginRow": { "type": "integer" },
              "endRow": { "type": "integer" },
              "tId": { "type": "integer" }
            },
            "required" : [ "beginRow", "endRow","tId" ]
          }
        },
        "version" : { "type": "string" }
      },
      "required": [
        "fileSections",
        "numReads",
        "basicData",
        "version"
      ]
    }

# Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
# Author: Derek Barnett

import PacBioBam
import config 

import os
import unittest

class EndToEndTest(unittest.TestCase):
        
    def runTest(self):
        
        testData = config.TestData()
        ex2BamFn = testData.directory + "/ex2.bam"
        generatedBamFn = testData.directory + "/generated.bam"
        
        # loop over original file, store names, write to generated file
        file = PacBioBam.BamFile(ex2BamFn)
        self.assertTrue(file)

        writer = PacBioBam.BamWriter(generatedBamFn, file.Header())
        self.assertTrue(writer)

        entireFile = PacBioBam.EntireFileQuery(file)
        self.assertTrue(entireFile)
         
        names_in = []
        for record in PacBioBam.Iterate(entireFile):
            names_in.append(record.FullName())
            writer.Write(record)

        writer.Close() # force flush & close
        
        # open generated file, read names
        file = PacBioBam.BamFile(generatedBamFn)
        self.assertTrue(file)
        
        entireFile = PacBioBam.EntireFileQuery(file)
        self.assertTrue(entireFile)
        
        names_out = []
        for record in PacBioBam.Iterate(entireFile):
            names_out.append(record.FullName())
        
        # ensure equal
        self.assertEqual(names_in, names_out)
        
        # clean up
        os.remove(generatedBamFn)

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
import unittest

class QualityValueTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def runTest(self):
        self.test_defaults()
        self.test_fromNumber()
        self.test_fromFastq()
        
    # ------------ TESTS --------------
    
    def test_defaults(self):
        value = PacBioBam.QualityValue()    
        self.assertEqual(0,  int(value))
        self.assertEqual('!', value.Fastq())
        
    def test_fromNumber(self):
        
        zero        = PacBioBam.QualityValue(0)
        thirtythree = PacBioBam.QualityValue(33)
        normal      = PacBioBam.QualityValue(42)
        maxQV       = PacBioBam.QualityValue(93)
        tooHigh     = PacBioBam.QualityValue(94)
        max8bit     = PacBioBam.QualityValue(126)
        
        self.assertEqual(0,  int(zero))
        self.assertEqual(33, int(thirtythree))
        self.assertEqual(42, int(normal))
        self.assertEqual(93, int(maxQV))
        self.assertEqual(93, int(tooHigh))
        self.assertEqual(93, int(max8bit))

        self.assertEqual('!', zero.Fastq())
        self.assertEqual('B', thirtythree.Fastq())
        self.assertEqual('K', normal.Fastq())
        self.assertEqual('~', maxQV.Fastq())
        self.assertEqual('~', tooHigh.Fastq())
        self.assertEqual('~', max8bit.Fastq())

    def test_fromFastq(self):
        
        zero        = PacBioBam.QualityValue.FromFastq('!')
        thirtythree = PacBioBam.QualityValue.FromFastq('B')
        normal      = PacBioBam.QualityValue.FromFastq('K')
        maxQV       = PacBioBam.QualityValue.FromFastq('~')
        
        self.assertEqual(0,  int(zero))
        self.assertEqual(33, int(thirtythree))
        self.assertEqual(42, int(normal))
        self.assertEqual(93, int(maxQV))
        
class QualityValuesTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def runTest(self):
        self.test_defaults()
        self.test_fromNumbers()
        self.test_fromFastq()
        
    # ------------ TESTS --------------
    
    def test_defaults(self):
        values = PacBioBam.QualityValues()    
        self.assertFalse(values.Fastq()) 
        
    def test_fromNumbers(self):
        
        fastqString = "~~~KKBB!!"
        values = [ 93, 93, 93, 42, 42, 33, 33, 0, 0 ]
        
        qvs = PacBioBam.QualityValues()
        for value in values:
            qvs.append(PacBioBam.QualityValue(value))
            
        self.assertEqual(fastqString, qvs.Fastq())

    def test_fromFastq(self):
        
        fastqString = "~~~KKBB!!"
        values = [ 93, 93, 93, 42, 42, 33, 33, 0, 0 ]
        
        qvs = PacBioBam.QualityValues.FromFastq(fastqString)
        
        self.assertEqual(len(fastqString), len(qvs))
        self.assertEqual(len(values),      len(qvs))
        
        for i, v in enumerate(values):
            self.assertEqual(v, int(qvs[i]))
    
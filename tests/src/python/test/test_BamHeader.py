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

class BamHeaderTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def runTest(self):
        self.test_defaultCtor()
        self.test_decode()
        self.test_encode()
        
    # ------------ TESTS --------------
    
    def test_defaultCtor(self):   
             
        header = PacBioBam.BamHeader()

        self.assertFalse(header.Version())
        self.assertFalse(header.SortOrder())
        self.assertEqual(0, len(header.ReadGroups()))
        self.assertEqual(0, len(header.Sequences()))
        self.assertEqual(0, len(header.Programs()))
        self.assertEqual(0, len(header.Comments()))
        
        with self.assertRaises(RuntimeError):
            pg = header.Program("foo")
            rg = header.ReadGroup("foo")
            sq = header.SequenceId("foo")
            sl = header.SequenceLength(42)
            sn = header.SequenceName(42)
        
        
    def test_decode(self):
        
        text = ("@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"
               "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
               "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
               "@RG\tID:rg1\tSM:control\n"
               "@RG\tID:rg2\tSM:condition1\n"
               "@RG\tID:rg3\tSM:condition1\n"
               "@PG\tID:_foo_\tPN:ide\n"
               "@CO\tipsum and so on\n"
               "@CO\tcitation needed\n")

        header = PacBioBam.BamHeader(text)

        self.assertEqual("1.1",       header.Version())
        self.assertEqual("queryname", header.SortOrder())
        self.assertEqual("3.0.1",     header.PacBioBamVersion())

        self.assertEqual(3, len(header.ReadGroups()))
        self.assertTrue(header.HasReadGroup("rg1"))
        self.assertTrue(header.HasReadGroup("rg2"))
        self.assertTrue(header.HasReadGroup("rg3"))
        self.assertEqual("control",    header.ReadGroup("rg1").Sample())
        self.assertEqual("condition1", header.ReadGroup("rg2").Sample())
        self.assertEqual("condition1", header.ReadGroup("rg3").Sample())

        self.assertEqual(2, len(header.Sequences()))
        self.assertTrue(header.HasSequence("chr1"))
        self.assertTrue(header.HasSequence("chr2"))
        self.assertEqual("chocobo", header.Sequence("chr1").Species())
        self.assertEqual("chocobo", header.Sequence("chr2").Species())
        self.assertEqual("2038",    header.Sequence("chr1").Length())
        self.assertEqual("3042",    header.Sequence("chr2").Length())

        self.assertEqual(1, len(header.Programs()))
        self.assertTrue(header.HasProgram("_foo_"))
        self.assertEqual("ide", header.Program("_foo_").Name())

        self.assertEqual(2, len(header.Comments()))
        self.assertEqual("ipsum and so on", header.Comments()[0])
        self.assertEqual("citation needed", header.Comments()[1])
        
    def test_encode(self):
        
        expectedText = ("@HD\tVN:1.1\tSO:queryname\tpb:3.0.1\n"
                        "@SQ\tSN:chr1\tLN:2038\tSP:chocobo\n"
                        "@SQ\tSN:chr2\tLN:3042\tSP:chocobo\n"
                        "@RG\tID:rg1\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:control\n"
                        "@RG\tID:rg2\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                        "@RG\tID:rg3\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tSM:condition1\n"
                        "@PG\tID:_foo_\tPN:ide\n"
                        "@CO\tipsum and so on\n"
                        "@CO\tcitation needed\n")

        rg1 = PacBioBam.ReadGroupInfo("rg1")
        rg1.Sample("control")
        rg2 = PacBioBam.ReadGroupInfo("rg2")
        rg2.Sample("condition1")
        rg3 = PacBioBam.ReadGroupInfo("rg3")
        rg3.Sample("condition1")

        seq1 = PacBioBam.SequenceInfo("chr1")
        seq1.Length("2038")
        seq1.Species("chocobo")
        seq2 = PacBioBam.SequenceInfo("chr2")
        seq2.Length("3042")
        seq2.Species("chocobo")

        prog1 = PacBioBam.ProgramInfo("_foo_")
        prog1.Name("ide")

        header = PacBioBam.BamHeader()
        header.Version("1.1")
        header.SortOrder("queryname")
        header.PacBioBamVersion("3.0.1")
        header.AddReadGroup(rg1)
        header.AddReadGroup(rg2)
        header.AddReadGroup(rg3)
        header.AddSequence(seq1)
        header.AddSequence(seq2)
        header.AddProgram(prog1)
        header.AddComment("ipsum and so on")
        header.AddComment("citation needed")

        self.assertEqual(expectedText, header.ToSam())
        
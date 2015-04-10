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

class CigarTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def runTest(self):
        self.test_typeToChar()
        self.test_charToType()
        self.test_setType()
        self.test_setChar()
        self.test_cigarOpCtors()
        self.test_fromEmptyString()
        self.test_fromString()
        self.test_toEmptyString()
        self.test_toString()
        
    # ------------ TESTS --------------
    
    def test_typeToChar(self):
        self.assertEqual('M', PacBioBam.CigarOperation.TypeToChar(PacBioBam.ALIGNMENT_MATCH))
        self.assertEqual('I', PacBioBam.CigarOperation.TypeToChar(PacBioBam.INSERTION))
        self.assertEqual('D', PacBioBam.CigarOperation.TypeToChar(PacBioBam.DELETION))
        self.assertEqual('N', PacBioBam.CigarOperation.TypeToChar(PacBioBam.REFERENCE_SKIP))
        self.assertEqual('S', PacBioBam.CigarOperation.TypeToChar(PacBioBam.SOFT_CLIP))
        self.assertEqual('H', PacBioBam.CigarOperation.TypeToChar(PacBioBam.HARD_CLIP))
        self.assertEqual('P', PacBioBam.CigarOperation.TypeToChar(PacBioBam.PADDING))
        self.assertEqual('=', PacBioBam.CigarOperation.TypeToChar(PacBioBam.SEQUENCE_MATCH))
        self.assertEqual('X', PacBioBam.CigarOperation.TypeToChar(PacBioBam.SEQUENCE_MISMATCH))
        
    def test_charToType(self):
        self.assertEqual(PacBioBam.ALIGNMENT_MATCH,   PacBioBam.CigarOperation.CharToType('M'))
        self.assertEqual(PacBioBam.INSERTION,         PacBioBam.CigarOperation.CharToType('I'))
        self.assertEqual(PacBioBam.DELETION,          PacBioBam.CigarOperation.CharToType('D'))
        self.assertEqual(PacBioBam.REFERENCE_SKIP,    PacBioBam.CigarOperation.CharToType('N'))
        self.assertEqual(PacBioBam.SOFT_CLIP,         PacBioBam.CigarOperation.CharToType('S'))
        self.assertEqual(PacBioBam.HARD_CLIP,         PacBioBam.CigarOperation.CharToType('H'))
        self.assertEqual(PacBioBam.PADDING,           PacBioBam.CigarOperation.CharToType('P'))
        self.assertEqual(PacBioBam.SEQUENCE_MATCH,    PacBioBam.CigarOperation.CharToType('='))
        self.assertEqual(PacBioBam.SEQUENCE_MISMATCH, PacBioBam.CigarOperation.CharToType('X'))
        
    def test_setType(self):
        m = PacBioBam.CigarOperation()
        i = PacBioBam.CigarOperation()
        d = PacBioBam.CigarOperation()
        n = PacBioBam.CigarOperation()
        s = PacBioBam.CigarOperation()
        h = PacBioBam.CigarOperation()
        p = PacBioBam.CigarOperation()
        e = PacBioBam.CigarOperation()
        x = PacBioBam.CigarOperation()
        
        m.Type(PacBioBam.ALIGNMENT_MATCH)
        i.Type(PacBioBam.INSERTION)
        d.Type(PacBioBam.DELETION)
        n.Type(PacBioBam.REFERENCE_SKIP)
        s.Type(PacBioBam.SOFT_CLIP)
        h.Type(PacBioBam.HARD_CLIP)
        p.Type(PacBioBam.PADDING)
        e.Type(PacBioBam.SEQUENCE_MATCH)
        x.Type(PacBioBam.SEQUENCE_MISMATCH)
        
        self.assertEqual('M', m.Char())
        self.assertEqual('I', i.Char())
        self.assertEqual('D', d.Char())
        self.assertEqual('N', n.Char())
        self.assertEqual('S', s.Char())
        self.assertEqual('H', h.Char())
        self.assertEqual('P', p.Char())
        self.assertEqual('=', e.Char())
        self.assertEqual('X', x.Char())
        
    def test_setChar(self):
        m = PacBioBam.CigarOperation()
        i = PacBioBam.CigarOperation()
        d = PacBioBam.CigarOperation()
        n = PacBioBam.CigarOperation()
        s = PacBioBam.CigarOperation()
        h = PacBioBam.CigarOperation()
        p = PacBioBam.CigarOperation()
        e = PacBioBam.CigarOperation()
        x = PacBioBam.CigarOperation()
        
        m.Char('M')
        i.Char('I')
        d.Char('D')
        n.Char('N')
        s.Char('S')
        h.Char('H')
        p.Char('P')
        e.Char('=')
        x.Char('X')
        
        self.assertEqual(PacBioBam.ALIGNMENT_MATCH,   m.Type())
        self.assertEqual(PacBioBam.INSERTION,         i.Type())
        self.assertEqual(PacBioBam.DELETION,          d.Type())
        self.assertEqual(PacBioBam.REFERENCE_SKIP,    n.Type())
        self.assertEqual(PacBioBam.SOFT_CLIP,         s.Type())
        self.assertEqual(PacBioBam.HARD_CLIP,         h.Type())
        self.assertEqual(PacBioBam.PADDING,           p.Type())
        self.assertEqual(PacBioBam.SEQUENCE_MATCH,    e.Type())
        self.assertEqual(PacBioBam.SEQUENCE_MISMATCH, x.Type())
        
    def test_cigarOpCtors(self):
        c1 = PacBioBam.CigarOperation('S', 10)
        c2 = PacBioBam.CigarOperation(PacBioBam.SOFT_CLIP, 10)
        
        self.assertEqual('S', c1.Char())
        self.assertEqual('S', c2.Char())
        self.assertEqual(PacBioBam.SOFT_CLIP, c1.Type())
        self.assertEqual(PacBioBam.SOFT_CLIP, c2.Type())
        self.assertEqual(10, c1.Length())
        self.assertEqual(10, c2.Length())
        
    def test_fromEmptyString(self):
        s = ""
        cigar = PacBioBam.Cigar(s)
        self.assertEqual(0, len(cigar))
        
    def test_fromString(self):
        singleCigarString = "100M"
        multiCigarString  = "100M2D34I6M"
        
        singleCigar = PacBioBam.Cigar(singleCigarString)
        multiCigar  = PacBioBam.Cigar(multiCigarString)
        
        self.assertEqual(1, len(singleCigar))
        c = singleCigar[0]
        self.assertEqual('M', c.Char())
        self.assertEqual(100, c.Length())

        self.assertEqual(4, len(multiCigar))
        op0 = multiCigar[0]
        op1 = multiCigar[1]
        op2 = multiCigar[2]
        op3 = multiCigar[3]

        self.assertEqual('M', op0.Char())
        self.assertEqual('D', op1.Char())
        self.assertEqual('I', op2.Char())
        self.assertEqual('M', op3.Char())
        self.assertEqual(100, op0.Length())
        self.assertEqual(2,   op1.Length())
        self.assertEqual(34,  op2.Length())
        self.assertEqual(6,   op3.Length())
        
    def test_toEmptyString(self):
        cigar = PacBioBam.Cigar()
        self.assertFalse(cigar.ToStdString())

    def test_toString(self):
        
        singleCigarString = "100M"
        multiCigarString  = "100M2D34I6M"
        
        singleCigar = PacBioBam.Cigar()
        singleCigar.append(PacBioBam.CigarOperation(PacBioBam.ALIGNMENT_MATCH, 100))
        
        multiCigar = PacBioBam.Cigar()
        multiCigar.append(PacBioBam.CigarOperation(PacBioBam.ALIGNMENT_MATCH, 100))
        multiCigar.append(PacBioBam.CigarOperation(PacBioBam.DELETION, 2))
        multiCigar.append(PacBioBam.CigarOperation(PacBioBam.INSERTION, 34))
        multiCigar.append(PacBioBam.CigarOperation(PacBioBam.ALIGNMENT_MATCH, 6))
        
        self.assertEqual(singleCigarString, singleCigar.ToStdString())
        self.assertEqual(multiCigarString,  multiCigar.ToStdString())
        
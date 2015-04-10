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

class IntervalsTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def runTest(self):
        self.test_unmappedPosition()
        self.test_ctors()
        self.test_equality()
        self.test_copy()
        self.test_modifiers()
        self.test_cover()
        self.test_intersect()
        self.test_validity()
        self.test_length()
        
    # ------------ TESTS --------------
    
    def test_unmappedPosition(self):
        self.assertEqual(-1, PacBioBam.UnmappedPosition)
    
    def test_ctors(self):
        empty  = PacBioBam.PositionInterval()
        single = PacBioBam.PositionInterval(4)
        normal = PacBioBam.PositionInterval(5, 8)
        
        self.assertEqual(0, empty.Start())
        self.assertEqual(0, empty.Stop())
        self.assertEqual(4, single.Start())
        self.assertEqual(5, single.Stop())
        self.assertEqual(5, normal.Start())
        self.assertEqual(8, normal.Stop())
        
    def test_equality(self):
        
        empty           = PacBioBam.PositionInterval()
        empty2          = PacBioBam.PositionInterval()
        singleton       = PacBioBam.PositionInterval(4)
        sameAsSingleton = PacBioBam.PositionInterval(4, 5)
        normal          = PacBioBam.PositionInterval(5, 8)
        sameAsNormal    = PacBioBam.PositionInterval(5, 8)
        different       = PacBioBam.PositionInterval(20, 40)
        
        # self-equality
        self.assertEqual(empty, empty)
        self.assertEqual(singleton, singleton)
        self.assertEqual(normal, normal)
        self.assertEqual(different, different)

        # same values
        self.assertEqual(empty, empty2)
        self.assertEqual(singleton, sameAsSingleton)
        self.assertEqual(normal, sameAsNormal)
        
        # different values
        self.assertNotEqual(empty, singleton)
        self.assertNotEqual(empty, normal)
        self.assertNotEqual(empty, different)
        self.assertNotEqual(singleton, normal)
        self.assertNotEqual(normal, different)
        
    def test_copy(self):
        interval1 = PacBioBam.PositionInterval(5,8)
        interval2 = PacBioBam.PositionInterval(interval1)
        interval3 = interval1
        
        self.assertEqual(interval1, interval1)
        self.assertEqual(interval1, interval2)
        self.assertEqual(interval1, interval3)
        
    def test_modifiers(self):
        
        interval1 = PacBioBam.PositionInterval(5,8)
        interval2 = PacBioBam.PositionInterval(interval1)
        interval2.Start(2)
        interval2.Stop(10)
        
        self.assertNotEqual(interval1, interval2)
        self.assertEqual(2,  interval2.Start())
        self.assertEqual(10, interval2.Stop())
        
    def test_cover(self):
        
        a = PacBioBam.PositionInterval(2,4)
        b = PacBioBam.PositionInterval(3,5)
        c = PacBioBam.PositionInterval(6,8)
        d = PacBioBam.PositionInterval(1,7)
        e = PacBioBam.PositionInterval(5,8)
        
        #   0123456789  
        # a   --
        # b    --
        # c       --
        # d  ------
        # e      ---

        # self-cover
        self.assertTrue(a.Covers(a)) 
        self.assertTrue(a.CoveredBy(a))

        # basic covers/covered
        self.assertTrue(b.CoveredBy(d))
        self.assertTrue(d.Covers(b))
        self.assertNotEqual(b, d)
        self.assertFalse(b.Covers(d))
 
        # completely disjoint
        self.assertFalse(b.Covers(c))
        self.assertFalse(c.Covers(b))
        self.assertFalse(b.CoveredBy(c))
        self.assertFalse(c.CoveredBy(b))
        
        # b.stop == e.start
        self.assertFalse(b.Covers(e))
        self.assertFalse(b.CoveredBy(e))

        # shared endpoint, start contained
        self.assertTrue(e.Covers(c))
        self.assertTrue(c.CoveredBy(e))
        
    def test_intersect(self):
        
        a = PacBioBam.PositionInterval(2,4)
        b = PacBioBam.PositionInterval(3,5)
        c = PacBioBam.PositionInterval(6,8)
        d = PacBioBam.PositionInterval(1,7)
        e = PacBioBam.PositionInterval(5,8)
        
        #   0123456789  
        # a   --
        # b    --
        # c       --
        # d  ------
        # e      ---
        
        # self-intersection
        self.assertTrue(a.Intersects(a))
        
        # intersection is commutative
        self.assertTrue(a.Intersects(b))
        self.assertTrue(b.Intersects(a))
        
        # covered implies intersection
        self.assertTrue(d.Covers(a))
        self.assertTrue(a.Intersects(d))
        self.assertTrue(d.Intersects(a))

        # c.start > b.stop (obvious disjoint)
        self.assertFalse(b.Intersects(c))
        
        # b.stop == e.start (intervals are right-open, so disjoint)
        self.assertFalse(b.Intersects(e))

    def test_validity(self):
        
        a = PacBioBam.PositionInterval()     # default ctor
        b = PacBioBam.PositionInterval(0,0)  # start == stop (zero)
        c = PacBioBam.PositionInterval(4,4)  # start == stop (nonzero)
        d = PacBioBam.PositionInterval(0,1)  # start < stop  (start is zero)
        e = PacBioBam.PositionInterval(4,5)  # start < stop  (start is nonzero)
        f = PacBioBam.PositionInterval(5,4)  # start > stop
        
        self.assertFalse(a.IsValid())
        self.assertFalse(b.IsValid())
        self.assertFalse(c.IsValid())
        self.assertTrue(d.IsValid())
        self.assertTrue(e.IsValid())
        self.assertFalse(f.IsValid())
        
    def test_length(self):
        
        a = PacBioBam.PositionInterval(2,4)
        b = PacBioBam.PositionInterval(3,5)
        c = PacBioBam.PositionInterval(6,8)
        d = PacBioBam.PositionInterval(1,7)
        e = PacBioBam.PositionInterval(5,8)
        
        self.assertEqual(2, a.Length())
        self.assertEqual(2, b.Length())
        self.assertEqual(2, c.Length())
        self.assertEqual(6, d.Length())
        self.assertEqual(3, e.Length())
        
class GenomicIntervalsTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def runTest(self):
        self.test_ctors()
        self.test_copy()
        self.test_modifiers()
        self.test_cover()
        self.test_validity()
    
    # ------------ TESTS --------------
    
    def test_ctors(self):
        
        empty  = PacBioBam.GenomicInterval()
        normal = PacBioBam.GenomicInterval(0, 100, 200)
        
        self.assertEqual(-1, empty.Id())
        self.assertEqual(0,  empty.Start())
        self.assertEqual(0,  empty.Stop())
        
        self.assertEqual(0,   normal.Id())
        self.assertEqual(100, normal.Start())
        self.assertEqual(200, normal.Stop())

        
    def test_copy(self):
        
        a = PacBioBam.GenomicInterval(1, 10, 20)
        b = PacBioBam.GenomicInterval(a)
        c = a
        
        self.assertEqual(a, a)
        self.assertEqual(a, b)
        self.assertEqual(a, c)
        
    def test_modifiers(self):
        
        a = PacBioBam.GenomicInterval(1, 10, 20)
        
        b = PacBioBam.GenomicInterval(a)
        b.Id(5).Start(2).Stop(10)
        
        c = PacBioBam.GenomicInterval(a)
        c.Interval(b.Interval())
        
        self.assertNotEqual(a, b)
        self.assertEqual(5,  b.Id())
        self.assertEqual(2,  b.Start())
        self.assertEqual(10, b.Stop())        
        self.assertEqual(a.Id(), c.Id())
        self.assertEqual(b.Interval(), c.Interval())
        
    def test_cover(self):
        
        a = PacBioBam.GenomicInterval(0,2,4)
        b = PacBioBam.GenomicInterval(0,3,5)
        c = PacBioBam.GenomicInterval(0,6,8)
        d = PacBioBam.GenomicInterval(0,1,7)
        e = PacBioBam.GenomicInterval(0,5,8)
        f = PacBioBam.GenomicInterval(1,3,5)  # same as b, different ref
        
        #   0123456789  
        # a   --
        # b    --
        # c       --
        # d  ------
        # e      ---
        
        # self-cover
        self.assertTrue(a.Covers(a))
        self.assertTrue(a.CoveredBy(a))
        
        # basic covers/covered
        self.assertTrue(b.CoveredBy(d))
        self.assertTrue(d.Covers(b))
        self.assertNotEqual(b, d)
        self.assertFalse(b.Covers(d))
        
        # same coords as b, but different ref
        self.assertFalse(f.CoveredBy(d))
        self.assertFalse(d.Covers(f))
        self.assertNotEqual(f, d)
        self.assertFalse(f.Covers(d))
        
        # obvious disjoint
        self.assertFalse(b.Covers(c))
        self.assertFalse(c.Covers(b))
        self.assertFalse(b.CoveredBy(c))
        self.assertFalse(c.CoveredBy(b))
        
        # b.stop == e.start (intervals are right-open, so disjoint)
        self.assertFalse(b.Covers(e))
        self.assertFalse(b.CoveredBy(e))
        
        # shared endpoint, start contained
        self.assertTrue(e.Covers(c))
        self.assertTrue(c.CoveredBy(e)) 
        
    def test_validity(self):
        
        a = PacBioBam.GenomicInterval()       # default
        b = PacBioBam.GenomicInterval(0,0,0)  # valid id, start == stop (zero)
        c = PacBioBam.GenomicInterval(0,4,4)  # valid id, start == stop (non-zero)
        d = PacBioBam.GenomicInterval(0,0,1)  # valid id, start <  stop (start == zero)     OK
        e = PacBioBam.GenomicInterval(0,4,5)  # valid id, start <  stop (start >  zero)     OK
        f = PacBioBam.GenomicInterval(0,5,4)  # valid id, start >  stop 
        g = PacBioBam.GenomicInterval(-1,0,0) # invalid id, start == stop (zero)
        h = PacBioBam.GenomicInterval(-1,4,4) # invalid id, start == stop (non-zero)
        i = PacBioBam.GenomicInterval(-1,0,1) # invalid id, start <  stop (start == zero)
        j = PacBioBam.GenomicInterval(-1,4,5) # invalid id, start <  stop (start >  zero)
        k = PacBioBam.GenomicInterval(-1,5,4) # invalid id, start >  stop 
             
        self.assertTrue(d.IsValid())
        self.assertTrue(e.IsValid())
        self.assertFalse(a.IsValid())
        self.assertFalse(b.IsValid())
        self.assertFalse(c.IsValid())
        self.assertFalse(f.IsValid())
        self.assertFalse(g.IsValid())
        self.assertFalse(h.IsValid())
        self.assertFalse(i.IsValid())
        self.assertFalse(j.IsValid())
        self.assertFalse(k.IsValid())

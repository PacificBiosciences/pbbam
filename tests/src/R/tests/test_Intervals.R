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
	
test_case("Intervals_UnmappedPosition", { 
	assertEqual(-1L, UnmappedPosition())
})

test_case("Intervals_Ctors", { 
	
	empty  <- PositionInterval()
	single <- PositionInterval(4)
	normal <- PositionInterval(5, 8)

    assertEqual(0L, empty$Start())
    assertEqual(0L, empty$Stop())
    assertEqual(4L, single$Start())
    assertEqual(5L, single$Stop())
    assertEqual(5L, normal$Start())
    assertEqual(8L, normal$Stop())
})

test_case("Intervals_Equality", { 
	
    empty           <- PositionInterval()
    empty2          <- PositionInterval()
    singleton       <- PositionInterval(4)
    sameAsSingleton <- PositionInterval(4, 5)
    normal          <- PositionInterval(5, 8)
    sameAsNormal    <- PositionInterval(5, 8)
    different       <- PositionInterval(20, 40)
    
    # self-equality
    assertEqual(empty,     empty)
    assertEqual(singleton, singleton)
    assertEqual(normal,    normal)
    assertEqual(different, different)

    # same values
	# TODO: fix this to work with == or *anything* cleaner
	assertTrue(empty$'__eq__'(empty2))
	assertTrue(singleton$'__eq__'(sameAsSingleton))
	assertTrue(normal$'__eq__'(sameAsNormal))
	
    # different values
    assertNotEqual(empty,     singleton)
    assertNotEqual(empty,     normal)
    assertNotEqual(empty,     different)
    assertNotEqual(singleton, normal)
    assertNotEqual(normal,    different)
})

test_case("Intervals_Copy", { 
	
    interval1 <- PositionInterval(5,8)
    interval2 <- PositionInterval(interval1)
    interval3 <- interval1

	# TODO: fix this to work with == or *anything* cleaner
	assertTrue(interval1$'__eq__'(interval1))
	assertTrue(interval1$'__eq__'(interval2))
	assertTrue(interval1$'__eq__'(interval3))
})

test_case("Intervals_Modifiers", { 
	
	interval1 <- PositionInterval(5,8)
    interval2 <- PositionInterval(interval1)
    interval2$Start(2)
    interval2$Stop(10)
        
    assertNotEqual(interval1, interval2)
    assertEqual(2L,  interval2$Start())
    assertEqual(10L, interval2$Stop())
})

test_case("Intervals_Cover", { 
	
    a <- PositionInterval(2,4)
    b <- PositionInterval(3,5)
    c <- PositionInterval(6,8)
    d <- PositionInterval(1,7)
    e <- PositionInterval(5,8)
    
    #   0123456789  
    # a   --
    # b    --
    # c       --
    # d  ------
    # e      ---

    # self-cover
    assertTrue(a$Covers(a)) 
    assertTrue(a$CoveredBy(a))

    # basic covers/covered
    assertTrue(b$CoveredBy(d))
    assertTrue(d$Covers(b))
    assertNotEqual(b, d)
    assertFalse(b$Covers(d))

    # completely disjoint
    assertFalse(b$Covers(c))
    assertFalse(c$Covers(b))
    assertFalse(b$CoveredBy(c))
    assertFalse(c$CoveredBy(b))
    
    # b.stop == e.start
    assertFalse(b$Covers(e))
    assertFalse(b$CoveredBy(e))

    # shared endpoint, start contained
    assertTrue(e$Covers(c))
    assertTrue(c$CoveredBy(e))
})

test_case("Intervals_Intersect", { 
	
    a <- PositionInterval(2,4)
    b <- PositionInterval(3,5)
    c <- PositionInterval(6,8)
    d <- PositionInterval(1,7)
    e <- PositionInterval(5,8)
    
    #   0123456789  
    # a   --
    # b    --
    # c       --
    # d  ------
    # e      ---
    
    # self-intersection
    assertTrue(a$Intersects(a))
    
    # intersection is commutative
    assertTrue(a$Intersects(b))
    assertTrue(b$Intersects(a))
    
    # covered implies intersection
    assertTrue(d$Covers(a))
    assertTrue(a$Intersects(d))
    assertTrue(d$Intersects(a))

    # c.start > b.stop (obvious disjoint)
    assertFalse(b$Intersects(c))
    
    # b.stop == e.start (intervals are right-open, so disjoint)
    assertFalse(b$Intersects(e))
})

test_case("Intervals_Validity", { 
	
    a <- PositionInterval()     # default ctor
    b <- PositionInterval(0,0)  # start == stop (zero)
    c <- PositionInterval(4,4)  # start == stop (nonzero)
    d <- PositionInterval(0,1)  # start < stop  (start is zero)
    e <- PositionInterval(4,5)  # start < stop  (start is nonzero)
    f <- PositionInterval(5,4)  # start > stop
    
    assertFalse(a$IsValid())
    assertFalse(b$IsValid())
    assertFalse(c$IsValid())
    assertTrue(d$IsValid())
    assertTrue(e$IsValid())
    assertFalse(f$IsValid())
})

test_case("Intervals_Length",{ 
	
    a <- PositionInterval(2,4)
    b <- PositionInterval(3,5)
    c <- PositionInterval(6,8)
    d <- PositionInterval(1,7)
    e <- PositionInterval(5,8)
    
    assertEqual(2L, a$Length())
    assertEqual(2L, b$Length())
    assertEqual(2L, c$Length())
    assertEqual(6L, d$Length())
    assertEqual(3L, e$Length()) 
})

test_case("GenomicIntervals_Ctors", { 
	
    empty  <- GenomicInterval()
    normal <- GenomicInterval("seq1", 100, 200)
    
    assertEqual("",  empty$Name())
    assertEqual(0L,  empty$Start())
    assertEqual(0L,  empty$Stop())
    
    assertEqual("seq1", normal$Name())
    assertEqual(100L,   normal$Start())
    assertEqual(200L,   normal$Stop())
})

test_case("GenomicIntervals_Copy", { 
	
    a <- GenomicInterval("seq1", 10, 20)
    b <- GenomicInterval(a)
    c <- a
    
	# TODO: fix this to work with == or *anything* cleaner
	assertTrue(a$'__eq__'(a))
	assertTrue(a$'__eq__'(b))
	assertTrue(a$'__eq__'(c))
})

test_case("GenomicIntervals_Modifiers", { 
	
    a <- GenomicInterval("seq1", 10, 20)
    
    b <- GenomicInterval(a)
    b$Name("seq5")
	b$Start(2)
	b$Stop(10)
    
    c <- GenomicInterval(a)
    c$Interval(b$Interval())
    
    assertNotEqual(a, b)
	
    assertEqual("seq5",  b$Name())
    assertEqual(2L,  b$Start())
    assertEqual(10L, b$Stop())        
	
    assertEqual(a$Name(), c$Name())
	
	# TODO: fix this to work with == or *anything* cleaner
	assertTrue(b$Interval()$'__eq__'(c$Interval()))
})

test_case("GenomicIntervals_Cover", { 
	
    a <- GenomicInterval("seq1",2,4)
    b <- GenomicInterval("seq1",3,5)
    c <- GenomicInterval("seq1",6,8)
    d <- GenomicInterval("seq1",1,7)
    e <- GenomicInterval("seq1",5,8)
    f <- GenomicInterval("seq2",3,5)  # same as b, different ref
    
    #   0123456789  
    # a   --
    # b    --
    # c       --
    # d  ------
    # e      ---
    
    # self-cover
    assertTrue(a$Covers(a))
    assertTrue(a$CoveredBy(a))
    
    # basic covers/covered
    assertTrue(b$CoveredBy(d))
    assertTrue(d$Covers(b))
    assertNotEqual(b, d)
    assertFalse(b$Covers(d))
    
    # same coords as b, but different ref
    assertFalse(f$CoveredBy(d))
    assertFalse(d$Covers(f))
    assertNotEqual(f, d)
    assertFalse(f$Covers(d))
    
    # obvious disjoint
    assertFalse(b$Covers(c))
    assertFalse(c$Covers(b))
    assertFalse(b$CoveredBy(c))
    assertFalse(c$CoveredBy(b))
    
    # b.stop == e.start (intervals are right-open, so disjoint)
    assertFalse(b$Covers(e))
    assertFalse(b$CoveredBy(e))
    
    # shared endpoint, start contained
    assertTrue(e$Covers(c))
    assertTrue(c$CoveredBy(e))
	
	# assertTrue(FALSE)
})

test_case("GenomicIntervals_Validity", { 
	
    a <- GenomicInterval()       # default
    b <- GenomicInterval("seq1",0,0)  # valid id, start == stop (zero)
    c <- GenomicInterval("seq1",4,4)  # valid id, start == stop (non-zero)
    d <- GenomicInterval("seq",0,1)  # valid id, start <  stop (start == zero)     OK
    e <- GenomicInterval("seq1",4,5)  # valid id, start <  stop (start >  zero)     OK
    f <- GenomicInterval("seq1",5,4)  # valid id, start >  stop 
    g <- GenomicInterval("",0,0) # invalid id, start == stop (zero)
    h <- GenomicInterval("",4,4) # invalid id, start == stop (non-zero)
    i <- GenomicInterval("",0,1) # invalid id, start <  stop (start == zero)
    j <- GenomicInterval("",4,5) # invalid id, start <  stop (start >  zero)
    k <- GenomicInterval("",5,4) # invalid id, start >  stop 
         
    assertTrue(d$IsValid())
    assertTrue(e$IsValid())
    assertFalse(a$IsValid())
    assertFalse(b$IsValid())
    assertFalse(c$IsValid())
    assertFalse(f$IsValid())
    assertFalse(g$IsValid())
    assertFalse(h$IsValid())
    assertFalse(i$IsValid())
    assertFalse(j$IsValid())
    assertFalse(k$IsValid())
})

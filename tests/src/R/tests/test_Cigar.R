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

test_case("Cigar_TypeToChar", {
    assertEqual('M', CigarOperation_TypeToChar('ALIGNMENT_MATCH'))
	assertEqual('I', CigarOperation_TypeToChar('INSERTION'))
	assertEqual('D', CigarOperation_TypeToChar('DELETION'))
	assertEqual('N', CigarOperation_TypeToChar('REFERENCE_SKIP'))
	assertEqual('S', CigarOperation_TypeToChar('SOFT_CLIP'))
	assertEqual('H', CigarOperation_TypeToChar('HARD_CLIP'))
	assertEqual('P', CigarOperation_TypeToChar('PADDING'))
	assertEqual('=', CigarOperation_TypeToChar('SEQUENCE_MATCH'))
	assertEqual('X', CigarOperation_TypeToChar('SEQUENCE_MISMATCH'))
})

test_case("Cigar_CharToType", {
	
    assertEqual('ALIGNMENT_MATCH',   CigarOperation_CharToType('M'))
    assertEqual('INSERTION',         CigarOperation_CharToType('I'))
    assertEqual('DELETION',          CigarOperation_CharToType('D'))
    assertEqual('REFERENCE_SKIP',    CigarOperation_CharToType('N'))
    assertEqual('SOFT_CLIP',         CigarOperation_CharToType('S'))
    assertEqual('HARD_CLIP',         CigarOperation_CharToType('H'))
    assertEqual('PADDING',           CigarOperation_CharToType('P'))
    assertEqual('SEQUENCE_MATCH',    CigarOperation_CharToType('='))
    assertEqual('SEQUENCE_MISMATCH', CigarOperation_CharToType('X'))
})

test_case("Cigar_SetType", {
	
    m = CigarOperation()
    i = CigarOperation()
    d = CigarOperation()
    n = CigarOperation()
    s = CigarOperation()
    h = CigarOperation()
    p = CigarOperation()
    e = CigarOperation()
    x = CigarOperation()
    
    m$Type('ALIGNMENT_MATCH')
    i$Type('INSERTION')
    d$Type('DELETION')
    n$Type('REFERENCE_SKIP')
    s$Type('SOFT_CLIP')
    h$Type('HARD_CLIP')
    p$Type('PADDING')
    e$Type('SEQUENCE_MATCH')
    x$Type('SEQUENCE_MISMATCH')
    
    assertEqual('M', m$Char())
    assertEqual('I', i$Char())
    assertEqual('D', d$Char())
    assertEqual('N', n$Char())
    assertEqual('S', s$Char())
    assertEqual('H', h$Char())
    assertEqual('P', p$Char())
    assertEqual('=', e$Char())
    assertEqual('X', x$Char())
})

test_case("Cigar_SetChar", {
	
    m = CigarOperation()
    i = CigarOperation()
    d = CigarOperation()
    n = CigarOperation()
    s = CigarOperation()
    h = CigarOperation()
    p = CigarOperation()
    e = CigarOperation()
    x = CigarOperation()

    m$Char('M')
    i$Char('I')
    d$Char('D')
    n$Char('N')
    s$Char('S')
    h$Char('H')
    p$Char('P')
    e$Char('=')
    x$Char('X')

    assertEqual('ALIGNMENT_MATCH',   m$Type())
    assertEqual('INSERTION',         i$Type())
    assertEqual('DELETION',          d$Type())
    assertEqual('REFERENCE_SKIP',    n$Type())
    assertEqual('SOFT_CLIP',         s$Type())
    assertEqual('HARD_CLIP',         h$Type())
    assertEqual('PADDING',           p$Type())
    assertEqual('SEQUENCE_MATCH',    e$Type())
    assertEqual('SEQUENCE_MISMATCH', x$Type())
})

test_case("Cigar_CigarOpCtors", {

	c1 <- CigarOperation("S", 10)
	c2 <- CigarOperation(CigarOperation_TypeToChar('SOFT_CLIP'), 10)

	assertEqual('S', c1$Char())
	assertEqual('S', c2$Char())
	assertEqual('SOFT_CLIP', c1$Type())
	assertEqual('SOFT_CLIP', c2$Type())
	assertEqual(10L, c1$Length())
	assertEqual(10L, c2$Length())
})

test_case("Cigar_FromEmptyString", {
	
	s <- ""
	cigar <- Cigar(s)
	assertEqual(0L, cigar$size())
})

test_case("Cigar_FromString", {
	
    singleCigarString <- "100M"
    multiCigarString  <- "100M2D34I6M"
    
    singleCigar <- Cigar(singleCigarString)
    multiCigar  <- Cigar(multiCigarString)
    
    assertEqual(1L, singleCigar$size())
	
	c <- singleCigar$front()
	assertEqual('M',               c$Char())
	assertEqual('ALIGNMENT_MATCH', c$Type())
	assertEqual(100L,              c$Length())

    assertEqual(4L, multiCigar$size())
	
	# haven't quite figured out [ ] accessors via SWIG, 
	# but this method does work w/ !ZERO!-based indices
    op0 <- multiCigar$'__getitem__'(0) 
    op1 <- multiCigar$'__getitem__'(1)
    op2 <- multiCigar$'__getitem__'(2)
    op3 <- multiCigar$'__getitem__'(3)
    
    assertEqual('M', op0$Char())
    assertEqual('D', op1$Char())
    assertEqual('I', op2$Char())
    assertEqual('M', op3$Char())
	assertEqual('ALIGNMENT_MATCH', op0$Type())
	assertEqual('DELETION',        op1$Type())
	assertEqual('INSERTION',       op2$Type())
	assertEqual('ALIGNMENT_MATCH', op3$Type())
    assertEqual(100L, op0$Length())
    assertEqual(2L,   op1$Length())
    assertEqual(34L,  op2$Length())
    assertEqual(6L,   op3$Length())
})

test_case("Cigar_ToEmptyString", {
	
	cigar <- Cigar()
	assertEqual(0L, nchar(cigar$ToStdString())) # empty string is 1
})

test_case("Cigar_ToString", {
	
    singleCigarString <- "100M"
    multiCigarString  <- "100M2D34I6M"
    
	singleCigar <- Cigar()
	singleCigar$push_back( CigarOperation(CigarOperation_TypeToChar('ALIGNMENT_MATCH'), 100) )

	multiCigar <- Cigar()
	multiCigar$push_back(CigarOperation('M', 100))
	multiCigar$push_back(CigarOperation('D', 2))
	multiCigar$push_back(CigarOperation('I', 34))
	multiCigar$push_back(CigarOperation('M', 6))

	assertEqual(singleCigarString, singleCigar$ToStdString())
	assertEqual(multiCigarString,  multiCigar$ToStdString())
})

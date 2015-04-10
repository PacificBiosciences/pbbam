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

test_case("QualityValue_Defaults", { 
	
	value <- QualityValue()
	assertEqual(0L, value$ToInt())
	assertEqual('!', value$Fastq())
})

test_case("QualityValue_FromNumber", { 
	
    zero        <- QualityValue(0)
    thirtythree <- QualityValue(33)
    normal      <- QualityValue(42)
    maxQV       <- QualityValue(93)
    tooHigh     <- QualityValue(94)
    max8bit     <- QualityValue(126)
    
    assertEqual(0L,  zero$ToInt())
    assertEqual(33L, thirtythree$ToInt())
    assertEqual(42L, normal$ToInt())
    assertEqual(93L, maxQV$ToInt())
    assertEqual(93L, tooHigh$ToInt())
    assertEqual(93L, max8bit$ToInt())

    assertEqual('!', zero$Fastq())
    assertEqual('B', thirtythree$Fastq())
    assertEqual('K', normal$Fastq())
    assertEqual('~', maxQV$Fastq())
    assertEqual('~', tooHigh$Fastq())
    assertEqual('~', max8bit$Fastq())
})

test_case("QualityValue_FromFastq", { 
	
    zero        <- QualityValue_FromFastq('!')
    thirtythree <- QualityValue_FromFastq('B')
    normal      <- QualityValue_FromFastq('K')
    maxQV       <- QualityValue_FromFastq('~')

    assertEqual(0L,  zero$ToInt())
    assertEqual(33L, thirtythree$ToInt())
    assertEqual(42L, normal$ToInt())
    assertEqual(93L, maxQV$ToInt())
})

test_case("QualityValues_Defaults", { 
    values <- QualityValues()   
	assertEqual(0L, nchar(values$Fastq()))
})

test_case("QualityValues_FromNumbers", { 
	
	fastqString <- "~~~KKBB!!"
	values <- c(93, 93, 93, 42, 42, 33, 33, 0, 0)
	
	qvs <- QualityValues()
	for (v in values) 
		qvs$push_back(QualityValue(v))
	
	assertEqual(fastqString, qvs$Fastq())
})

test_case("QualityValues_FromFastq", { 
	
	fastqString <- "~~~KKBB!!"
	values <- c(93L, 93L, 93L, 42L, 42L, 33L, 33L, 0L, 0L)
	
	qvs <- QualityValues(fastqString)
	assertEqual(nchar(fastqString), qvs$size())
	assertEqual(length(values),     qvs$size())
	
	numValues <- length(values)
	for ( i in 1:numValues ) {
		qv <- qvs$'__getitem__'(i-1)
		assertEqual(values[i], qv$ToInt())
	}
})

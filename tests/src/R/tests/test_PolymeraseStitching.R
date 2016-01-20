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
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES LOSS OF
# USE, DATA, OR PROFITS OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
# Author: Derek Barnett

compareContainers <- function(c1, c2) {
    
	assertEqual(length(c1), length(c2))
    
	numElements <- length(c1)
	for (i in 1:numElements)
		assertEqual(c1[i], c2[i])
}

compareFrames <- function(f1, f2) {
    
    d1 <- f1$Data()
    d2 <- f2$Data()
    compareContainers(d1, d2)
}

compareRecords <- function(b1, b2) {

    assertTrue(b1$HasDeletionQV())
    assertTrue(b1$HasDeletionTag())
    assertTrue(b1$HasInsertionQV())
    assertTrue(b1$HasMergeQV())
    assertTrue(b1$HasSubstitutionQV())
    assertTrue(b1$HasSubstitutionTag())
    assertTrue(b1$HasLabelQV())
    assertTrue(b1$HasAltLabelQV())
    assertTrue(b1$HasAltLabelTag())
    assertTrue(b1$HasPkmean())
    assertTrue(b1$HasPkmid())
    assertTrue(b1$HasPulseCall())
    assertTrue(b1$HasIPD())
    assertTrue(b1$HasPulseWidth())
    assertTrue(b1$HasPrePulseFrames())
    assertTrue(b1$HasPulseCallWidth())
    assertTrue(b1$HasPulseMergeQV())
    
    assertTrue(b2$HasDeletionQV())
    assertTrue(b2$HasDeletionTag())
    assertTrue(b2$HasInsertionQV())
    assertTrue(b2$HasMergeQV())
    assertTrue(b2$HasSubstitutionQV())
    assertTrue(b2$HasSubstitutionTag())
    assertTrue(b2$HasLabelQV())
    assertTrue(b2$HasAltLabelQV())
    assertTrue(b2$HasAltLabelTag())
    assertTrue(b2$HasPkmean())
    assertTrue(b2$HasPkmid())
    assertTrue(b2$HasPulseCall())
    assertTrue(b2$HasIPD())
    assertTrue(b2$HasPulseWidth())
    assertTrue(b2$HasPrePulseFrames())
    assertTrue(b2$HasPulseCallWidth())
    assertTrue(b2$HasPulseMergeQV())
 
    assertEqual(b1$FullName(),        b2$FullName())
    assertEqual(b1$HoleNumber(),      b2$HoleNumber())
    assertEqual(b1$NumPasses(),       b2$NumPasses())
    assertEqual(b1$Sequence(),        b2$Sequence())
    assertEqual(b1$DeletionTag(),     b2$DeletionTag())
    assertEqual(b1$SubstitutionTag(), b2$SubstitutionTag())
    assertEqual(b1$AltLabelTag(),     b2$AltLabelTag())
    assertEqual(b1$PulseCall(),       b2$PulseCall())
    
    # compareContainers(b1$Pkmean(), b2$Pkmean())
    # compareContainers(b1$Pkmid(), b2$Pkmid())
    #
    # compareFrames(b1$IPD(),             b2$IPD())
    # compareFrames(b1$PulseWidth(),      b2$PulseWidth())
    # compareFrames(b1$PrePulseFrames(),  b2$PrePulseFrames())
    # compareFrames(b1$PulseCallWidth(),  b2$PulseCallWidth())

    assertEqual(b1$ReadGroup()$Id(), b2$ReadGroup()$Id())
    
    assertEqual(b1$Qualities()$Fastq(),       b2$Qualities()$Fastq())
    assertEqual(b1$DeletionQV()$Fastq(),      b2$DeletionQV()$Fastq())
    assertEqual(b1$InsertionQV()$Fastq(),     b2$InsertionQV()$Fastq())
    assertEqual(b1$MergeQV()$Fastq(),         b2$MergeQV()$Fastq())
    assertEqual(b1$SubstitutionQV()$Fastq(),  b2$SubstitutionQV()$Fastq())
    assertEqual(b1$LabelQV()$Fastq(),         b2$LabelQV()$Fastq())
    assertEqual(b1$AltLabelQV()$Fastq(),      b2$AltLabelQV()$Fastq())
    assertEqual(b1$PulseMergeQV()$Fastq(),    b2$PulseMergeQV()$Fastq())
    
    return
}

getVirtualRecord <- function(fn1, fn2) {
    
    result <- tryCatch(
        {
            vpr <- VirtualPolymeraseReader(fn1, fn2)
            
            assertTrue(vpr$HasNext())
            
            virtualRecord <- vpr$Next()
            
            assertFalse(vpr$HasNext())
            
            return(virtualRecord)
        },
        error = function(e) {
            print(paste('e:',e))
            assertTrue(FALSE) # should not throw
            return
        }    
    )
    return(result)
}

getPolymeraseRecord <- function(fn) {
    
    result <- tryCatch(
        {
            ds <- DataSet(fn)
            entireFile <- EntireFileQuery(ds)
            
            polyIter <- entireFile$begin()
            polyEnd <- entireFile$end()
            
            assertTrue(polyIter$'__ne__'(polyEnd))
            
            polyRecord <- polyIter$value()
            polyIter$incr()
            
            assertTrue(polyIter$'__eq__'(polyEnd))
            
            return(polyRecord)
        },
        error = function(e) {
            print(paste('e:',e))
            assertTrue(FALSE) # should not throw
            return
        }    
    )
    return(result)
}

test_case("PolymeraseStitching_VirtualRegions", {
	
	subreadsFn <- paste(test_data_path, "polymerase/internal.subreads.bam", sep="/")
	scrapsFn   <- paste(test_data_path, "polymerase/internal.scraps.bam", sep="/")
    virtualRecord <- getVirtualRecord(subreadsFn, scrapsFn)
    
    # -- ADAPTER -- #
    
    adapter <- virtualRecord$VirtualRegionsTable('ADAPTER')
    assertEqual(7L, adapter$size())
    
    region <- adapter$'__getitem__'(0)
    assertEqual(3047L, region$beginPos)
    assertEqual(3095L, region$endPos)
    
    region <- adapter$'__getitem__'(1)
    assertEqual(3650L, region$beginPos)
    assertEqual(3700L, region$endPos)
    
    region <- adapter$'__getitem__'(2)
    assertEqual(4289L, region$beginPos)
    assertEqual(4335L, region$endPos)
    
    region <- adapter$'__getitem__'(3)
    assertEqual(4888L, region$beginPos)
    assertEqual(4939L, region$endPos)
    
    region <- adapter$'__getitem__'(4)
    assertEqual(5498L, region$beginPos)
    assertEqual(5546L, region$endPos)
    
    region <- adapter$'__getitem__'(5)
    assertEqual(6116L, region$beginPos)
    assertEqual(6173L, region$endPos)
    
    region <- adapter$'__getitem__'(6)
    assertEqual(6740L, region$beginPos)
    assertEqual(6790L, region$endPos)

    # -- BARCODE -- #

    barcode = virtualRecord$VirtualRegionsTable('BARCODE')
    assertEqual(14L, barcode$size())

    region <- barcode$'__getitem__'(0)
    assertEqual(3025L, region$beginPos)
    assertEqual(3047L, region$endPos)
    
    region <- barcode$'__getitem__'(1)
    assertEqual(3095L, region$beginPos)
    assertEqual(3116L, region$endPos)
    
    region <- barcode$'__getitem__'(2)
    assertEqual(3628L, region$beginPos)
    assertEqual(3650L, region$endPos)
    
    region <- barcode$'__getitem__'(3)
    assertEqual(3700L, region$beginPos)
    assertEqual(3722L, region$endPos)
    
    region <- barcode$'__getitem__'(4)
    assertEqual(4267L, region$beginPos)
    assertEqual(4289L, region$endPos)
    
    region <- barcode$'__getitem__'(5)
    assertEqual(4335L, region$beginPos)
    assertEqual(4356L, region$endPos)
    
    region <- barcode$'__getitem__'(6)
    assertEqual(4864L, region$beginPos)
    assertEqual(4888L, region$endPos)

    region <- barcode$'__getitem__'(7)
    assertEqual(4939L, region$beginPos)
    assertEqual(4960L, region$endPos)
    
    region <- barcode$'__getitem__'(8)
    assertEqual(5477L, region$beginPos)
    assertEqual(5498L, region$endPos)
    
    region <- barcode$'__getitem__'(9)
    assertEqual(5546L, region$beginPos)
    assertEqual(5571L, region$endPos)
    
    region <- barcode$'__getitem__'(10)
    assertEqual(6087L, region$beginPos)
    assertEqual(6116L, region$endPos)
    
    region <- barcode$'__getitem__'(11)
    assertEqual(6173L, region$beginPos)
    assertEqual(6199L, region$endPos)
    
    region <- barcode$'__getitem__'(12)
    assertEqual(6719L, region$beginPos)
    assertEqual(6740L, region$endPos)
    
    region <- barcode$'__getitem__'(13)
    assertEqual(6790L, region$beginPos)
    assertEqual(6812L, region$endPos)

    # -- LQREGION -- #

    lqregion = virtualRecord$VirtualRegionsTable('LQREGION')
    assertEqual(2L, lqregion$size())
    
    region <- lqregion$'__getitem__'(0)
    assertEqual(0L, region$beginPos)
    assertEqual(2659L, region$endPos)
    
    region <- lqregion$'__getitem__'(1)
    assertEqual(7034L, region$beginPos)
    assertEqual(7035L, region$endPos)
    
    # -- HQREGION -- #

    hqregion = virtualRecord$VirtualRegionsTable('HQREGION')
    assertEqual(1L, hqregion$size())
    
    region <- hqregion$'__getitem__'(0)
    assertEqual(2659L, region$beginPos)
    assertEqual(7034L, region$endPos)
})

test_case("PolymeraseStitching_InternalSubreadsToOriginal", {
  
    # stitch virtual polymerase record
    subreadsFn <- paste(test_data_path, "polymerase/internal.subreads.bam", sep="/")
    scrapsFn   <- paste(test_data_path, "polymerase/internal.scraps.bam", sep="/")
    virtualRecord <- getVirtualRecord(subreadsFn, scrapsFn)

    # fetch original polymerase record
    polyFn <- paste(test_data_path, "polymerase/internal.polymerase.bam", sep="/")
    polyRecord <- getPolymeraseRecord(polyFn)      

    # check
    compareRecords(polyRecord, virtualRecord)
})

test_case("PolymeraseStitching_InternalHQToOriginal", {
  
    # stitch virtual polymerase record
    hqRegionFn <- paste(test_data_path, "polymerase/internal.hqregions.bam", sep="/")
    lqRegionFn <- paste(test_data_path, "polymerase/internal.lqregions.bam", sep="/")
    virtualRecord <- getVirtualRecord(hqRegionFn, lqRegionFn)
    
    # fetch original polymerase record
    polyFn <- paste(test_data_path, "polymerase/internal.polymerase.bam", sep="/")
    polyRecord <- getPolymeraseRecord(polyFn)      

    # check
    compareRecords(polyRecord, virtualRecord)
})

test_case("PolymeraseStitching_ProductionSubreadsToOriginal", {
  
    # stitch virtual polymerase record
    subreadsFn <- paste(test_data_path, "polymerase/production.subreads.bam", sep="/")
    scrapsFn   <- paste(test_data_path, "polymerase/production.scraps.bam", sep="/")
    virtualRecord <- getVirtualRecord(subreadsFn, scrapsFn)
    
    # fetch original polymerase record
    polyFn <- paste(test_data_path, "polymerase/production.polymerase.bam", sep="/")
    polyRecord <- getPolymeraseRecord(polyFn)  
    
    # compare
    assertEqual(polyRecord$FullName(),        virtualRecord$FullName())
    assertEqual(polyRecord$HoleNumber(),      virtualRecord$HoleNumber())
    assertEqual(polyRecord$NumPasses(),       virtualRecord$NumPasses())
    assertEqual(polyRecord$Sequence(),        virtualRecord$Sequence())
    assertEqual(polyRecord$DeletionTag(),     virtualRecord$DeletionTag())
    assertEqual(polyRecord$SubstitutionTag(), virtualRecord$SubstitutionTag())
    
    compareFrames(polyRecord$IPD(),                virtualRecord$IPDV1Frames())
    assertEqual(polyRecord$ReadGroup()$Id(),       virtualRecord$ReadGroup()$Id())
    
    tolerance = 1e-5
    assertTrue( abs(polyRecord$ReadAccuracy()$ToFloat() - virtualRecord$ReadAccuracy()$ToFloat()) <= tolerance )
    # assertEqual(polyRecord$ReadAccuracy()$ToFloat(), virtualRecord$ReadAccuracy()$ToFloat())

    assertEqual(polyRecord$Qualities()$Fastq(),       virtualRecord$Qualities()$Fastq())
    assertEqual(polyRecord$DeletionQV()$Fastq(),      virtualRecord$DeletionQV()$Fastq())
    assertEqual(polyRecord$InsertionQV()$Fastq(),     virtualRecord$InsertionQV()$Fastq())
    assertEqual(polyRecord$MergeQV()$Fastq(),         virtualRecord$MergeQV()$Fastq())
    assertEqual(polyRecord$SubstitutionQV()$Fastq(),  virtualRecord$SubstitutionQV()$Fastq())
})

test_case("PolymeraseStitching_ProductionHQToOriginal", {
  
    # stitch virtual polymerase record
    hqRegionFn <- paste(test_data_path, "polymerase/production_hq.hqregion.bam", sep="/")
    lqRegionFn <- paste(test_data_path, "polymerase/production_hq.scraps.bam", sep="/")
    virtualRecord <- getVirtualRecord(hqRegionFn, lqRegionFn)
  
    # fetch original polymerase record
    polyFn <- paste(test_data_path, "polymerase/production.polymerase.bam", sep="/")
    polyRecord <- getPolymeraseRecord(polyFn)
  
    # compare
    assertEqual(polyRecord$FullName(),        virtualRecord$FullName())
    assertEqual(polyRecord$HoleNumber(),      virtualRecord$HoleNumber())
    assertEqual(polyRecord$NumPasses(),       virtualRecord$NumPasses())
    assertEqual(polyRecord$Sequence(),        virtualRecord$Sequence())
    assertEqual(polyRecord$DeletionTag(),     virtualRecord$DeletionTag())
    assertEqual(polyRecord$SubstitutionTag(), virtualRecord$SubstitutionTag())

    compareFrames(polyRecord$IPD(),                virtualRecord$IPDV1Frames())
    assertEqual(polyRecord$ReadGroup()$Id(),       virtualRecord$ReadGroup()$Id())
    
    tolerance = 1e-5
    assertTrue( abs(polyRecord$ReadAccuracy()$ToFloat() - virtualRecord$ReadAccuracy()$ToFloat()) <= tolerance )
    # assertEqual(polyRecord$ReadAccuracy()$ToInt(), virtualRecord$ReadAccuracy()$ToInt())
    
    assertEqual(polyRecord$Qualities()$Fastq(),       virtualRecord$Qualities()$Fastq())
    assertEqual(polyRecord$DeletionQV()$Fastq(),      virtualRecord$DeletionQV()$Fastq())
    assertEqual(polyRecord$InsertionQV()$Fastq(),     virtualRecord$InsertionQV()$Fastq())
    assertEqual(polyRecord$MergeQV()$Fastq(),         virtualRecord$MergeQV()$Fastq())
    assertEqual(polyRecord$SubstitutionQV()$Fastq(),  virtualRecord$SubstitutionQV()$Fastq())

    assertTrue(polyRecord$HasDeletionQV())
    assertTrue(polyRecord$HasDeletionTag())
    assertTrue(polyRecord$HasInsertionQV())
    assertTrue(polyRecord$HasMergeQV())
    assertTrue(polyRecord$HasSubstitutionQV())
    assertTrue(polyRecord$HasSubstitutionTag())
    assertTrue(polyRecord$HasIPD())
    assertFalse(polyRecord$HasLabelQV())
    assertFalse(polyRecord$HasAltLabelQV())
    assertFalse(polyRecord$HasAltLabelTag())
    assertFalse(polyRecord$HasPkmean())
    assertFalse(polyRecord$HasPkmid())
    assertFalse(polyRecord$HasPulseCall())
    assertFalse(polyRecord$HasPulseWidth())
    assertFalse(polyRecord$HasPrePulseFrames())
    assertFalse(polyRecord$HasPulseCallWidth())
    assertFalse(polyRecord$HasPulseCall())

    assertTrue(virtualRecord$HasDeletionQV())
    assertTrue(virtualRecord$HasDeletionTag())
    assertTrue(virtualRecord$HasInsertionQV())
    assertTrue(virtualRecord$HasMergeQV())
    assertTrue(virtualRecord$HasSubstitutionQV())
    assertTrue(virtualRecord$HasSubstitutionTag())
    assertTrue(virtualRecord$HasIPD())
    assertFalse(virtualRecord$HasLabelQV())
    assertFalse(virtualRecord$HasAltLabelQV())
    assertFalse(virtualRecord$HasAltLabelTag())
    assertFalse(virtualRecord$HasPkmean())
    assertFalse(virtualRecord$HasPkmid())
    assertFalse(virtualRecord$HasPulseCall())
    assertFalse(virtualRecord$HasPulseWidth())
    assertFalse(virtualRecord$HasPrePulseFrames())
    assertFalse(virtualRecord$HasPulseCallWidth()) 
    assertFalse(virtualRecord$HasPulseCall())
})
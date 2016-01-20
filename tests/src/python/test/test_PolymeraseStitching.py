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

class PolymeraseStitchingTest(unittest.TestCase):
    
    # ------------ SETUP --------------
    
    def setUp(self):
        self.data = config.TestData()
    
    def runTest(self):
        self.test_virtualRegions()
        self.test_internalSubreadsToOriginal()
        self.test_internalHqToOriginal()
        self.test_productionSubreadsToOriginal()
        self.test_productionHqToOriginal()
        
    # ------------ TESTS --------------
    
    def test_virtualRegions(self):

        subreadBam = self.data.directory + "/polymerase/internal.subreads.bam"
        scrapsBam  = self.data.directory + "/polymerase/internal.scraps.bam"
        vpr = PacBioBam.VirtualPolymeraseReader(subreadBam, scrapsBam)
    
        virtualRecord = vpr.Next()
        
        # NOTE: this method is disabled 
        #
        # Any attempt to retrive this value resulted in several 
        #   "swig/python detected a memory leak of type 'unknown', no destructor found."
        # errors (& an empty dictionary result). The same info is available via the 
        # VirtualRegionsTable(regionType) method, though a bit clunkier if you just want 
        # to iterate. But access to region info for specific types are available & correct, 
        # so I'm just going to leave this one out for now. - DB
        #
        # regionMap = virtualRecord.VirtualRegionsMap();
    
        # ADAPTER
        adapter = virtualRecord.VirtualRegionsTable(PacBioBam.VirtualRegionType_ADAPTER)
        self.assertEqual(7, len(adapter))
        self.assertEqual(3047, adapter[0].beginPos);
        self.assertEqual(3095, adapter[0].endPos);
        self.assertEqual(3650, adapter[1].beginPos);
        self.assertEqual(3700, adapter[1].endPos);
        self.assertEqual(4289, adapter[2].beginPos);
        self.assertEqual(4335, adapter[2].endPos);
        self.assertEqual(4888, adapter[3].beginPos);
        self.assertEqual(4939, adapter[3].endPos);
        self.assertEqual(5498, adapter[4].beginPos);
        self.assertEqual(5546, adapter[4].endPos);
        self.assertEqual(6116, adapter[5].beginPos);
        self.assertEqual(6173, adapter[5].endPos);
        self.assertEqual(6740, adapter[6].beginPos);
        self.assertEqual(6790, adapter[6].endPos);
    
        # BARCODE
        barcode = virtualRecord.VirtualRegionsTable(PacBioBam.VirtualRegionType_BARCODE)
        self.assertEqual(14, len(barcode))
        self.assertEqual(3025, barcode[0].beginPos);
        self.assertEqual(3047, barcode[0].endPos);
        self.assertEqual(3095, barcode[1].beginPos);
        self.assertEqual(3116, barcode[1].endPos);
        self.assertEqual(3628, barcode[2].beginPos);
        self.assertEqual(3650, barcode[2].endPos);
        self.assertEqual(3700, barcode[3].beginPos);
        self.assertEqual(3722, barcode[3].endPos);
        self.assertEqual(4267, barcode[4].beginPos);
        self.assertEqual(4289, barcode[4].endPos);
        self.assertEqual(4335, barcode[5].beginPos);
        self.assertEqual(4356, barcode[5].endPos);
        self.assertEqual(4864, barcode[6].beginPos);
        self.assertEqual(4888, barcode[6].endPos);
        self.assertEqual(4939, barcode[7].beginPos);
        self.assertEqual(4960, barcode[7].endPos);
        self.assertEqual(5477, barcode[8].beginPos);
        self.assertEqual(5498, barcode[8].endPos);
        self.assertEqual(5546, barcode[9].beginPos);
        self.assertEqual(5571, barcode[9].endPos);
        self.assertEqual(6087, barcode[10].beginPos);
        self.assertEqual(6116, barcode[10].endPos);
        self.assertEqual(6173, barcode[11].beginPos);
        self.assertEqual(6199, barcode[11].endPos);
        self.assertEqual(6719, barcode[12].beginPos);
        self.assertEqual(6740, barcode[12].endPos);
        self.assertEqual(6790, barcode[13].beginPos);
        self.assertEqual(6812, barcode[13].endPos);
    
        # HQREGION
        hqregion = virtualRecord.VirtualRegionsTable(PacBioBam.VirtualRegionType_HQREGION)
        self.assertEqual(1, len(hqregion))
        
        self.assertEqual(2659, hqregion[0].beginPos);
        self.assertEqual(7034, hqregion[0].endPos);
    
        # LQREGION
        lqregion = virtualRecord.VirtualRegionsTable(PacBioBam.VirtualRegionType_LQREGION)
        self.assertEqual(2, len(lqregion))
        
        self.assertEqual(0,    lqregion[0].beginPos);
        self.assertEqual(2659, lqregion[0].endPos);
        self.assertEqual(7034, lqregion[1].beginPos);
        self.assertEqual(7035, lqregion[1].endPos);
    
        # SUBREAD
        subread = virtualRecord.VirtualRegionsTable(PacBioBam.VirtualRegionType_SUBREAD)
        self.assertEqual(8, len(subread)) 
    
    def test_internalSubreadsToOriginal(self):
        
        # stitch virtual polymerase record
        subreadsBam = self.data.directory + "/polymerase/internal.subreads.bam"
        scrapsBam   = self.data.directory + "/polymerase/internal.scraps.bam"
        vpr = PacBioBam.VirtualPolymeraseReader(subreadsBam, scrapsBam)

        self.assertTrue(vpr.HasNext())
        virtualRecord = vpr.Next()
        self.assertFalse(vpr.HasNext())

        # fetch original polymerase record
        polyBam   = PacBioBam.DataSet(self.data.directory + "/polymerase/internal.polymerase.bam")
        polyQuery = PacBioBam.EntireFileQuery(polyBam)

        polyIter = polyQuery.begin()
        polyEnd  = polyQuery.end()

        self.assertTrue(polyIter != polyEnd)
        polyRecord = polyIter.value()
        polyIter.incr()
        self.assertTrue(polyIter == polyEnd)

        # compare
        self.compare(polyRecord, virtualRecord)

    def test_internalHqToOriginal(self):
        
        # stitch virtual polymerase record
        hqRegionsBam = self.data.directory + "/polymerase/internal.hqregions.bam"
        lqRegionsBam = self.data.directory + "/polymerase/internal.lqregions.bam"
        vpr = PacBioBam.VirtualPolymeraseReader(hqRegionsBam, lqRegionsBam)

        self.assertTrue(vpr.HasNext())
        virtualRecord = vpr.Next()
        self.assertFalse(vpr.HasNext())

        # fetch original polymerase record
        polyBam   = PacBioBam.DataSet(self.data.directory + "/polymerase/internal.polymerase.bam")
        polyQuery = PacBioBam.EntireFileQuery(polyBam)

        polyIter = polyQuery.begin()
        polyEnd  = polyQuery.end()

        self.assertTrue(polyIter != polyEnd)
        polyRecord = polyIter.value()
        polyIter.incr()
        self.assertTrue(polyIter == polyEnd)
       
        # # compare
        self.compare(polyRecord, virtualRecord)

    def test_productionSubreadsToOriginal(self):
        
        # stitch virtual polymerase record
        subreadsBam = self.data.directory + "/polymerase/production.subreads.bam"
        scrapsBam   = self.data.directory + "/polymerase/production.scraps.bam"
        vpr = PacBioBam.VirtualPolymeraseReader(subreadsBam, scrapsBam)

        self.assertTrue(vpr.HasNext())
        virtualRecord = vpr.Next()
        self.assertFalse(vpr.HasNext())

        # fetch original polymerase record
        polyBam   = PacBioBam.DataSet(self.data.directory + "/polymerase/production.polymerase.bam")
        polyQuery = PacBioBam.EntireFileQuery(polyBam)

        polyIter = polyQuery.begin()
        polyEnd  = polyQuery.end()

        self.assertTrue(polyIter != polyEnd)
        polyRecord = polyIter.value()
        polyIter.incr()
        self.assertTrue(polyIter == polyEnd)
        
        # compare
        self.assertEqual(polyRecord.FullName(),        virtualRecord.FullName());
        self.assertEqual(polyRecord.HoleNumber(),      virtualRecord.HoleNumber());
        self.assertEqual(polyRecord.NumPasses(),       virtualRecord.NumPasses());
        self.assertEqual(polyRecord.Sequence(),        virtualRecord.Sequence());
        self.assertEqual(polyRecord.DeletionTag(),     virtualRecord.DeletionTag());
        self.assertEqual(polyRecord.SubstitutionTag(), virtualRecord.SubstitutionTag());
        self.assertEqual(polyRecord.IPD(),             virtualRecord.IPDV1Frames());
        self.assertEqual(polyRecord.ReadGroup(),       virtualRecord.ReadGroup());
        
        self.assertAlmostEqual(float(polyRecord.ReadAccuracy()), float(virtualRecord.ReadAccuracy()));
        
        self.assertEqual(polyRecord.Qualities().Fastq(),       virtualRecord.Qualities().Fastq());
        self.assertEqual(polyRecord.DeletionQV().Fastq(),      virtualRecord.DeletionQV().Fastq());
        self.assertEqual(polyRecord.InsertionQV().Fastq(),     virtualRecord.InsertionQV().Fastq());
        self.assertEqual(polyRecord.MergeQV().Fastq(),         virtualRecord.MergeQV().Fastq());
        self.assertEqual(polyRecord.SubstitutionQV().Fastq(),  virtualRecord.SubstitutionQV().Fastq());

    def test_productionHqToOriginal(self):
        
        # stitch virtual polymerase record
        hqRegionsBam = self.data.directory + "/polymerase/production_hq.hqregion.bam"
        lqRegionsBam = self.data.directory + "/polymerase/production_hq.scraps.bam"
        vpr = PacBioBam.VirtualPolymeraseReader(hqRegionsBam, lqRegionsBam)

        self.assertTrue(vpr.HasNext())
        virtualRecord = vpr.Next()
        self.assertFalse(vpr.HasNext())

        # fetch original polymerase record
        polyBam   = PacBioBam.DataSet(self.data.directory + "/polymerase/production.polymerase.bam")
        polyQuery = PacBioBam.EntireFileQuery(polyBam)

        polyIter = polyQuery.begin()
        polyEnd  = polyQuery.end()

        self.assertTrue(polyIter != polyEnd)
        polyRecord = polyIter.value()
        polyIter.incr()
        self.assertTrue(polyIter == polyEnd)
        
        # compare        
        self.assertFalse(polyRecord.HasPulseCall());
        self.assertFalse(virtualRecord.HasPulseCall());
        
        self.assertEqual(polyRecord.FullName(),        virtualRecord.FullName());
        self.assertEqual(polyRecord.HoleNumber(),      virtualRecord.HoleNumber());
        self.assertEqual(polyRecord.NumPasses(),       virtualRecord.NumPasses());
        self.assertEqual(polyRecord.Sequence(),        virtualRecord.Sequence());
        self.assertEqual(polyRecord.DeletionTag(),     virtualRecord.DeletionTag());
        self.assertEqual(polyRecord.SubstitutionTag(), virtualRecord.SubstitutionTag());
        self.assertEqual(polyRecord.IPD(),             virtualRecord.IPDV1Frames());
        self.assertEqual(polyRecord.ReadGroup(),       virtualRecord.ReadGroup());
        
        self.assertAlmostEqual(float(polyRecord.ReadAccuracy()), float(virtualRecord.ReadAccuracy()));
        
        self.assertEqual(polyRecord.Qualities().Fastq(),       virtualRecord.Qualities().Fastq());
        self.assertEqual(polyRecord.DeletionQV().Fastq(),      virtualRecord.DeletionQV().Fastq());
        self.assertEqual(polyRecord.InsertionQV().Fastq(),     virtualRecord.InsertionQV().Fastq());
        self.assertEqual(polyRecord.MergeQV().Fastq(),         virtualRecord.MergeQV().Fastq());
        self.assertEqual(polyRecord.SubstitutionQV().Fastq(),  virtualRecord.SubstitutionQV().Fastq());
        
        self.assertTrue(polyRecord.HasDeletionQV());
        self.assertTrue(polyRecord.HasDeletionTag());
        self.assertTrue(polyRecord.HasInsertionQV());
        self.assertTrue(polyRecord.HasMergeQV());
        self.assertTrue(polyRecord.HasSubstitutionQV());
        self.assertTrue(polyRecord.HasSubstitutionTag());
        self.assertTrue(polyRecord.HasIPD());
        self.assertFalse(polyRecord.HasLabelQV());
        self.assertFalse(polyRecord.HasAltLabelQV());
        self.assertFalse(polyRecord.HasAltLabelTag());
        self.assertFalse(polyRecord.HasPkmean());
        self.assertFalse(polyRecord.HasPkmid());
        self.assertFalse(polyRecord.HasPulseCall());
        self.assertFalse(polyRecord.HasPulseWidth());
        self.assertFalse(polyRecord.HasPrePulseFrames());
        self.assertFalse(polyRecord.HasPulseCallWidth());
        
        self.assertTrue(virtualRecord.HasDeletionQV());
        self.assertTrue(virtualRecord.HasDeletionTag());
        self.assertTrue(virtualRecord.HasInsertionQV());
        self.assertTrue(virtualRecord.HasMergeQV());
        self.assertTrue(virtualRecord.HasSubstitutionQV());
        self.assertTrue(virtualRecord.HasSubstitutionTag());
        self.assertTrue(virtualRecord.HasIPD());
        self.assertFalse(virtualRecord.HasLabelQV());
        self.assertFalse(virtualRecord.HasAltLabelQV());
        self.assertFalse(virtualRecord.HasAltLabelTag());
        self.assertFalse(virtualRecord.HasPkmean());
        self.assertFalse(virtualRecord.HasPkmid());
        self.assertFalse(virtualRecord.HasPulseCall());
        self.assertFalse(virtualRecord.HasPulseWidth());
        self.assertFalse(virtualRecord.HasPrePulseFrames());
        self.assertFalse(virtualRecord.HasPulseCallWidth());   
    
    # ------------ HELPERS --------------
    
    def compare(self, b1, b2):
    
        self.assertTrue(b1.HasDeletionQV());
        self.assertTrue(b1.HasDeletionTag());
        self.assertTrue(b1.HasInsertionQV());
        self.assertTrue(b1.HasMergeQV());
        self.assertTrue(b1.HasSubstitutionQV());
        self.assertTrue(b1.HasSubstitutionTag());
        self.assertTrue(b1.HasLabelQV());
        self.assertTrue(b1.HasAltLabelQV());
        self.assertTrue(b1.HasAltLabelTag());
        self.assertTrue(b1.HasPkmean());
        self.assertTrue(b1.HasPkmid());
        self.assertTrue(b1.HasPulseCall());
        self.assertTrue(b1.HasIPD());
        self.assertTrue(b1.HasPulseWidth());
        self.assertTrue(b1.HasPrePulseFrames());
        self.assertTrue(b1.HasPulseCallWidth());
        self.assertTrue(b1.HasPulseMergeQV());

        self.assertTrue(b2.HasDeletionQV());
        self.assertTrue(b2.HasDeletionTag());
        self.assertTrue(b2.HasInsertionQV());
        self.assertTrue(b2.HasMergeQV());
        self.assertTrue(b2.HasSubstitutionQV());
        self.assertTrue(b2.HasSubstitutionTag());
        self.assertTrue(b2.HasLabelQV());
        self.assertTrue(b2.HasAltLabelQV());
        self.assertTrue(b2.HasAltLabelTag());
        self.assertTrue(b2.HasPkmean());
        self.assertTrue(b2.HasPkmid());
        self.assertTrue(b2.HasPulseCall());
        self.assertTrue(b2.HasIPD());
        self.assertTrue(b2.HasPulseWidth());
        self.assertTrue(b2.HasPrePulseFrames());
        self.assertTrue(b2.HasPulseCallWidth());
        self.assertTrue(b2.HasPulseMergeQV());
    
        self.assertEqual(b1.FullName(),        b2.FullName());
        self.assertEqual(b1.HoleNumber(),      b2.HoleNumber());
        self.assertEqual(b1.NumPasses(),       b2.NumPasses());
        self.assertEqual(b1.Sequence(),        b2.Sequence());
        self.assertEqual(b1.DeletionTag(),     b2.DeletionTag());
        self.assertEqual(b1.SubstitutionTag(), b2.SubstitutionTag());
        self.assertEqual(b1.AltLabelTag(),     b2.AltLabelTag());
        self.assertEqual(b1.Pkmean(),          b2.Pkmean());
        self.assertEqual(b1.Pkmid(),           b2.Pkmid());
        self.assertEqual(b1.PulseCall(),       b2.PulseCall());
        self.assertEqual(b1.IPD(),             b2.IPD());
        self.assertEqual(b1.PulseWidth(),      b2.PulseWidth());
        self.assertEqual(b1.PrePulseFrames(),  b2.PrePulseFrames());
        self.assertEqual(b1.PulseCallWidth(),  b2.PulseCallWidth());
        self.assertEqual(b1.ReadGroup(),       b2.ReadGroup());
        
        self.assertEqual(b1.Qualities().Fastq(),       b2.Qualities().Fastq());
        self.assertEqual(b1.DeletionQV().Fastq(),      b2.DeletionQV().Fastq());
        self.assertEqual(b1.InsertionQV().Fastq(),     b2.InsertionQV().Fastq());
        self.assertEqual(b1.MergeQV().Fastq(),         b2.MergeQV().Fastq());
        self.assertEqual(b1.SubstitutionQV().Fastq(),  b2.SubstitutionQV().Fastq());
        self.assertEqual(b1.PulseMergeQV().Fastq(),    b2.PulseMergeQV().Fastq());
        self.assertEqual(b1.LabelQV().Fastq(),         b2.LabelQV().Fastq());
        self.assertEqual(b1.AltLabelQV().Fastq(),      b2.AltLabelQV().Fastq());
        
        

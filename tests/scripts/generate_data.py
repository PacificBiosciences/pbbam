#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals

import os, shutil, sys
from io import StringIO

# FASTA generation
fastaSeq_1 = """TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
AACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGCAGGCG"""

fastaSeq_2 = """TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
AACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAAC"""

fastaSeq_3 = """TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
ACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCT"""

# FASTQ generation

fastqSeq_1   = """TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAACTCCGCCGGCGCAGGCG"""
fastqQuals_1 = """[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["""

fastqSeq_2   = """TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCTGTGCTGAGGAGAACGCAAC"""
fastqQuals_2 = """[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["""

fastqSeq_3   = """TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACACCCTAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCTAACCCCTAACCCTAACCCT"""
fastqQuals_3 = """]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]"""


# file creation decorator
def fileMaker(func):
    def inner(*args, **kwargs):
        print(" - Creating file: %s..." % args[1], end='')
        sys.stdout.flush()
        retval = func(*args)
        print("done.")
        sys.stdout.flush()
        return retval
    return inner

# symlink creation decorator
def fileLinker(func):
    def inner(*args, **kwargs):
        print(" - Creating symlink: %s..." % args[1], end='')
        sys.stdout.flush()
        retval = func(*args)
        print("done.")
        sys.stdout.flush()
        return retval
    return inner

# return a copy of original, minues any lines that contain an entry in blacklist
def trimXmlElements(original, blacklist):
    out = StringIO()
    for line in original.splitlines():
        if all(x not in line for x in blacklist):
            out.write(line + '\n')
    result = out.getvalue()
    out.close()
    return result

class TestDataGenerator:

    def __init__(self, source, dest):

        # source/destination directories
        self.testDataDir      = source
        self.generatedDataDir = dest

        # generated output files/symlinks & 'maker' functions
        self.outputFiles = {
            'truncated.bam' : self.makeTruncatedBam,
            'chunking_emptyfilters.subreadset.xml'   : self.makeChunkingXml,
            'chunking_missingfilters.subreadset.xml' : self.makeChunkingXml,
            'normal.fa' : self.makeNormalFasta,
            'normal.fq' : self.makeNormalFastq
        }
        self.outputSymlinks = {
            'aligned.bam'      : self.makeAlignedBamCopy,
            'aligned.bam.bai'  : self.makeAlignedBamCopy,
            'aligned.bam.pbi'  : self.makeAlignedBamCopy,
            'aligned2.bam'     : self.makeAlignedBamCopy,
            'aligned2.bam.bai' : self.makeAlignedBamCopy,
            'aligned2.bam.pbi' : self.makeAlignedBamCopy,
            'm150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam'     : self.makeChunkingSymlink,
            'm150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam.pbi' : self.makeChunkingSymlink,
            'm150404_101626_42267_c100807920800000001823174110291514_s1_p0.2.subreads.bam'     : self.makeChunkingSymlink,
            'm150404_101626_42267_c100807920800000001823174110291514_s1_p0.2.subreads.bam.pbi' : self.makeChunkingSymlink,
            'm150404_101626_42267_c100807920800000001823174110291514_s1_p0.3.subreads.bam'     : self.makeChunkingSymlink,
            'm150404_101626_42267_c100807920800000001823174110291514_s1_p0.3.subreads.bam.pbi' : self.makeChunkingSymlink,
            'missing_pbi.bam' : self.makeMissingPbiBam,
        }

    def editChunkingXml(self, outputFn, removeFiltersNode):
        inputXmlFn  = os.path.join(self.testDataDir,'chunking','chunking.subreadset.xml')
        outputXmlFn = os.path.join(self.generatedDataDir,outputFn)

        blacklist = ['pbds:Filter>', 'pbbase:Properties>', '<pbbase:Property']
        if removeFiltersNode:
            blacklist.append('pbds:Filters>')

        inputXml = ''
        with open(inputXmlFn, 'r') as xml_infile:
            inputXml = xml_infile.read()
        outputXml = trimXmlElements(inputXml, blacklist)
        with open(outputXmlFn, 'w') as xml_outfile:
            xml_outfile.write(outputXml)

    @fileLinker
    def makeAlignedBamCopy(self, outputFn):
        source = os.path.join(self.testDataDir,outputFn)
        dest   = os.path.join(self.generatedDataDir, outputFn)
        os.symlink(source, dest)

    @fileLinker
    def makeChunkingSymlink(self, outputFn):
        source = os.path.join(self.testDataDir,'chunking', outputFn)
        dest   = os.path.join(self.generatedDataDir, outputFn)
        os.symlink(source, dest)
  
    @fileLinker
    def makeMissingPbiBam(self, outputFn):
        source = os.path.join(self.testDataDir, 'phi29.bam')
        dest   = os.path.join(self.generatedDataDir, outputFn)
        os.symlink(source, dest)

    @fileMaker
    def makeChunkingXml(self, outputFn):
        if outputFn == 'chunking_emptyfilters.subreadset.xml':
            removeFiltersNode = False
        else:
            removeFiltersNode = True
        self.editChunkingXml(outputFn, removeFiltersNode)

    @fileMaker
    def makeNormalFasta(self, outputFn):
        content = ">1\n" + fastaSeq_1 + "\n>2\n" + fastaSeq_2 + "\n>3\n" + fastaSeq_3
        dest = os.path.join(self.generatedDataDir, outputFn)
        with open(outputFn, 'w') as fasta_out:
            fasta_out.write(content)

    @fileMaker
    def makeNormalFastq(self, outputFn):
        content = ("@1\n" + fastqSeq_1 + "\n+\n" + fastqQuals_1 + "\n" +
                   "@2\n" + fastqSeq_2 + "\n+\n" + fastqQuals_2 + "\n" +
                   "@3\n" + fastqSeq_3 + "\n+\n" + fastqQuals_3 + "\n")
        dest = os.path.join(self.generatedDataDir, outputFn)
        with open(outputFn, 'w') as fastq_out:
            fastq_out.write(content)

    @fileMaker
    def makeTruncatedBam(self, outputFn):
        source = os.path.join(self.testDataDir, 'phi29.bam')
        dest   = os.path.join(self.generatedDataDir, outputFn)
        shutil.copyfile(source, dest)
        with open(dest, 'r+b') as in_file:
            in_file.truncate(200)

    # main entry point
    def generate(self):

        # skip file if it exists
        os.chdir(self.generatedDataDir)
        filenames = list(self.outputFiles.keys())
        for file in filenames:
            if os.path.exists(file) :
                del self.outputFiles[file]

        # skip symlink if it exists
        symlinks = list(self.outputSymlinks.keys())
        for link in symlinks:
            if os.path.lexists(link):
                del self.outputSymlinks[link]

        # only print message & run makers, if any files/symlinks to be created
        # else silent success
        if self.outputFiles or self.outputSymlinks:
            print('Generating test data in %s ' % self.generatedDataDir)
            for file, func in self.outputFiles.items():
                func(file)
            for link, func in self.outputSymlinks.items():
                func(link)

# script entry point
if __name__ == '__main__':
    g = TestDataGenerator(sys.argv[1], sys.argv[2])
    g.generate()

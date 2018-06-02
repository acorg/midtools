from os.path import join
from math import log10
from collections import Counter

from dark.fasta import FastaReads

from data.component import connectedComponentsByOffset
from data.utils import baseCountsToStr


class ReadAnalysis(object):
    """
    Hold an entire read analysis.
    """
    def __init__(self, fastaFile, outputDir, genomeLength, alignedReads,
                 readCountAtOffset, baseCountAtOffset, readsAtOffset,
                 significantOffsets, threshold):
        self.fastaFile = fastaFile
        self.outputDir = outputDir
        self.genomeLength = genomeLength
        self.alignedReads = alignedReads
        self.readCountAtOffset = readCountAtOffset
        self.baseCountAtOffset = baseCountAtOffset
        self.readsAtOffset = readsAtOffset
        self.significantOffsets = significantOffsets
        self.threshold = threshold

        significantReads = set(read for read in alignedReads
                               if read.significantOffsets)

        # Find all connected components.
        self.components = components = []
        for count, component in enumerate(
                connectedComponentsByOffset(significantReads, threshold),
                start=1):
            component.saveFasta(outputDir, count)
            components.append(component)
            with open(join(outputDir, 'component-%d.txt' % count), 'w') as fp:
                component.summarize(fp, count)

        # Sanity check: The significantReads set should be be empty
        # following the above processing.
        assert len(significantReads) == 0

    def saveSignificantOffsets(self):
        with open(join(self.outputDir, 'significant-offsets.txt'), 'w') as fp:
            for offset in self.significantOffsets:
                print(offset, file=fp)

    def saveBaseFrequencies(self):
        genomeLengthWidth = int(log10(self.genomeLength)) + 1
        nucleotides = set('ACGT')

        with open(join(self.outputDir, 'base-frequencies.txt'), 'w') as fp:
            for offset in range(self.genomeLength):
                counts = Counter()
                for read in self.readsAtOffset[offset]:
                    base = read.base(offset)
                    if base in nucleotides:
                        counts[base] += 1
                print('Location %*d: base counts %s' %
                      (genomeLengthWidth, offset + 1, baseCountsToStr(counts)),
                      file=fp)

    def saveReducedFasta(self):
        """
        Write out FASTA that contains reads with bases just at the
        significant offsets.
        """
        allGaps = '-' * len(self.significantOffsets)

        def unwanted(read):
            return (None if read.sequence == allGaps else read)

        FastaReads(self.fastaFile).filter(
            keepSites=self.significantOffsets).filter(
                modifier=unwanted).save(join(self.outputDir, 'reduced.fasta'))

    def summarize(self):
        with open(join(self.outputDir, 'component-summary.txt'), 'w') as fp:

            print('Read %d aligned reads of length %d. '
                  'Found %d significant locations.' %
                  (len(self.alignedReads), self.genomeLength,
                   len(self.significantOffsets)), file=fp)

            print('Reads were assigned to %d components:' %
                  len(self.components), file=fp)

            totalReads = 0
            for i, component in enumerate(self.components, start=1):
                componentCount = len(component)
                totalReads += componentCount
                ccCounts = sorted(
                    map(len, (cc.reads
                              for cc in component.consistentComponents)),
                    reverse=True)
                print('component %d (%d reads), offsets %d to %d: %r' % (
                    i, componentCount,
                    min(component.offsets), max(component.offsets),
                    ccCounts), end='', file=fp)
                if len(ccCounts) > 1:
                    print(' ratio %.2f' % (ccCounts[0] / ccCounts[1]), file=fp)
                else:
                    print()

            # Sanity check.
            assert totalReads == len(self.alignedReads)

            print('In total, %d reads were assigned to components.' %
                  totalReads, file=fp)

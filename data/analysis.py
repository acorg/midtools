from os.path import join
from os import unlink
from math import log10
from collections import Counter
from pathlib import Path  # This is Python 3 only.
from itertools import chain
from pprint import pprint

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.reads import ReadsInRAM, Read, Reads

from data.component import connectedComponentsByOffset
from data.utils import baseCountsToStr, commonest


class ReadAnalysis(object):
    """
    Hold an entire read analysis.
    """
    def __init__(self, fastaFile, outputDir, genomeLength, alignedReads,
                 readCountAtOffset, baseCountAtOffset, readsAtOffset,
                 significantOffsets, threshold, verbose, referenceGenomeFiles):
        self.fastaFile = fastaFile
        self.outputDir = outputDir
        self.genomeLength = genomeLength
        self.alignedReads = alignedReads
        self.readCountAtOffset = readCountAtOffset
        self.baseCountAtOffset = baseCountAtOffset
        self.readsAtOffset = readsAtOffset
        self.significantOffsets = significantOffsets
        self.threshold = threshold
        self.verbose = verbose
        self._removePreExistingOutputFiles()
        self.referenceGenomes = self._readReferenceGenomes(
            referenceGenomeFiles)
        self.components = self._findConnectedComponents()

    def _readReferenceGenomes(self, referenceGenomeFiles):
        """
        Read reference genomes from files and check their lengths match the
        aligned reads.
        """
        result = ReadsInRAM()
        if referenceGenomeFiles:
            for filename in referenceGenomeFiles:
                for read in FastaReads(filename):
                    if len(read) == self.genomeLength:
                        result.add(read)
                    else:
                        print('WARNING: Reference genome %s in file %s has '
                              'length %d but the provided aligned reads were '
                              'aligned against something of length %d. This '
                              'reference will be ignored.' % (
                                  read.id, filename, len(read),
                                  self.genomeLength))
            if self.verbose:
                print('Read %d reference genomes' % len(result))
        return result

    def _findConnectedComponents(self):
        """
        Find all connected components.
        """
        significantReads = set(read for read in self.alignedReads
                               if read.significantOffsets)
        components = []
        for count, component in enumerate(
                connectedComponentsByOffset(significantReads, self.threshold),
                start=1):
            components.append(component)

        # Sanity check: The significantReads set should be be empty
        # following the above processing.
        assert len(significantReads) == 0
        return components

    def _removePreExistingOutputFiles(self):
        # Remove all pre-existing files from the output directory.
        # This prevents us from doing a run that results in (say) 6
        # component files and then later doing a run that results in
        # only 5 components and erroneously thinking that
        # component-6-2.fasta etc. are from the most recent run.
        dir = self.outputDir

        paths = list(map(str, chain(
            Path(dir).glob('component-*.fasta'),
            Path(dir).glob('component-*.txt'),
            Path(dir).glob('significant-offsets.txt'))))

        if paths:
            if self.verbose:
                print('Removing %d pre-existing output file%s from %s '
                      'directory.' %
                      (len(paths), '' if len(paths) == 1 else 's', dir))
            list(map(unlink, paths))

    def saveComponentFasta(self):
        for count, component in enumerate(self.components):
            component.saveFasta(self.outputDir, count)

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

    def saveClosestReferenceConsensuses(self):
        """
        Calculate and save the best consensus to each reference genome.
        """
        def ccMatchCount(cc, reference, drawFp, drawMessage):
            """
            Count the matches between a consistent component and a reference
            genome.

            @param cc: A C{ConsistentComponent} instance.
            @param reference: A C{Read} instance.
            @param drawFp: A file pointer to write information about draws (if
                any) to.
            @param drawMessage: A C{str} message to write to C{drawFp}. If the
                string contains '%(baseCounts)s' that will be replaced by a
                string representation of the base counts (in C{counts})
                obtained from C{baseCountsToStr}. If not, the base count info
                will be printed after the message.
            @return: The C{int} count of bases that match the reference
                for the offsets covered by the consistent component.
            """
            refSeq = reference.sequence
            nucleotides = cc.nucleotides
            count = 0
            for offset in nucleotides:
                message = drawMessage + (
                    ' offset %d: base counts' % offset) + ' %(baseCounts)s.'
                componentBase = commonest(nucleotides[offset], drawFp, message)
                count += int(componentBase == refSeq[offset])
            return count

        def bestConsistentComponent(component, reference, fp):
            """
            Find the best consistent component in the given
            C{ComponentByOffsets} instance

            @param component: A C{ComponentByOffsets} instance.
            @param reference: A C{Read} instance.
            @param fp: A file pointer to write information to.
            @return: The C{int} offset of the best consistent component.
            """
            bestScore = -1
            bestIndex = None
            offsetCount = len(component.offsets)
            for index, cc in enumerate(component.consistentComponents):
                # To compute how good each consistent component of a
                # ComponentByOffsets instance is, it's not enough to just
                # count the matches (or the fraction of matches) in the
                # consistent component because those components can be very
                # small (e.g., with just one read that may only cover one
                # offset) and with a perfect (1.0) internal match fraction.
                #
                # So we compute a score that is the product of 1) the
                # fraction of matches within the consistent component and
                # 2) the fraction of the ComponentByOffsets offsets that
                # are covered by the consistent component. A consistent
                # component that agrees perfectly with the reference at all
                # its covered offsets and which covers all the offset in
                # the ComponentByOffsets will have a score of 1.0
                matchCount = ccMatchCount(
                    cc, reference, fp,
                    '    Consistent component %d base draw' % index)
                print('  Consistent component %d (%d reads) has %d exact '
                      'matches with the reference, out of the %d offsets it '
                      'covers (%.2f%%).'
                      % (index, len(cc.reads), matchCount, len(cc.nucleotides),
                         matchCount / len(cc.nucleotides) * 100.0),
                      file=fp)
                score = matchCount / offsetCount
                if score == bestScore:
                    print('    WARNING: score %.2f draw with consistent '
                          'component %d' % (score, index), file=fp)
                elif score > bestScore:
                    bestScore = score
                    bestIndex = index

            print('  The best consistent component is index %d.' % bestIndex,
                  file=fp)

            return bestIndex

        baseCountAtOffset = self.baseCountAtOffset

        for reference in self.referenceGenomes:

            fields = reference.id.split(maxsplit=1)
            if len(fields) == 1:
                shortReferenceId = reference.id
                referenceIdRest = ''
            else:
                shortReferenceId = fields[0]
                referenceIdRest = ' ' + fields[1]

            infoFile = join(self.outputDir,
                            '%s-consensus.txt' % shortReferenceId)
            with open(infoFile, 'w') as infoFp:
                offsetsDone = set()
                consensus = [None] * self.genomeLength
                bestCcIndices = []
                for count, component in enumerate(self.components, start=1):
                    print('\nExamimining component %d with %d offsets: %s' %
                          (count, len(component.offsets),
                           ', '.join(map(str, sorted(component.offsets)))),
                          file=infoFp)
                    bestCcIndex = bestConsistentComponent(component, reference,
                                                          infoFp)
                    bestCcIndices.append(bestCcIndex)
                    bestCc = component.consistentComponents[bestCcIndex]
                    print('  Adding best nucleotides to consensus:', file=infoFp)
                    for offset in bestCc.nucleotides:
                        assert consensus[offset] is None
                        base = commonest(
                            bestCc.nucleotides[offset], infoFp,
                            ('    WARNING: base count draw at offset %d ' %
                             offset) + ' %(baseCounts)s.')
                        referenceBase = reference.sequence[offset]
                        if base == referenceBase:
                            mismatch = ''
                        else:
                            consensusBase = commonest(
                                baseCountAtOffset[offset], infoFp,
                                ('    WARNING: consensus base count draw at '
                                 'offset %d ' % offset) + ' %(baseCounts)s.')
                            mismatch = (
                                ' (mismatch: reference has %s, all-read '
                                'consensus has %s)' % (referenceBase,
                                                       consensusBase))

                        print('    Offset %d: %s from nucleotides %s%s' %
                              (offset, base,
                               baseCountsToStr(bestCc.nucleotides[offset]),
                               mismatch), file=infoFp)

                        consensus[offset] = base
                        offsetsDone.add(offset)

                # Fill in (from the overall read consensus) the offsets
                # that were not significant in any connected component.
                consensusOffsets = set(range(self.genomeLength)) - offsetsDone
                print('\nAdding bases from %d non-connected-component '
                      'consensus offsets:' % len(consensusOffsets),
                      file=infoFp)
                for offset in consensusOffsets:
                    assert consensus[offset] is None
                    base = commonest(
                        baseCountAtOffset[offset], infoFp,
                        ('    WARNING: consensus base count draw at offset %d '
                         % offset) + ' %(baseCounts)s.')
                    print('  Offset %d: %s from nucleotides %s' %
                          (offset, base,
                           baseCountsToStr(baseCountAtOffset[offset])),
                          file=infoFp, end='')
                    referenceBase = reference.sequence[offset]
                    if base == referenceBase:
                        print(file=infoFp)
                    else:
                        print(' (mismatch: reference has %s)' %
                              referenceBase, file=infoFp)
                    consensus[offset] = base

                filename = join(self.outputDir,
                                '%s-consensus.fasta' % shortReferenceId)
                consensusId = (
                    '%s-consensus best-consistent-components:%s%s' %
                    (shortReferenceId, ','.join(map(str, bestCcIndices)),
                     referenceIdRest))

                consensusRead = Read(consensusId, ''.join(consensus))

                match = compareDNAReads(reference, consensusRead)

                identity = (
                    match['match']['identicalMatchCount'] +
                    match['match']['ambiguousMatchCount']) / self.genomeLength

                print('\nIdentity of reference with consensus: identity %.2f'
                      % (identity * 100.0), file=infoFp)
                print('\nMatch info:', file=infoFp)
                pprint(match, indent=4, stream=infoFp)
                Reads([reference, consensusRead]).save(filename)

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

    def saveComponentConsensuses(self):
        for count, component in enumerate(self.components, start=1):
            component.saveConsensuses(self.outputDir, count)

    def summarize(self):
        with open(join(self.outputDir, 'component-summary.txt'), 'w') as fp:

            print('Read %d aligned reads of length %d. '
                  'Found %d significant locations.' %
                  (len(self.alignedReads), self.genomeLength,
                   len(self.significantOffsets)), file=fp)

            print('Reads were assigned to %d connected components:' %
                  len(self.components), file=fp)

            totalReads = 0
            for count, component in enumerate(self.components, start=1):

                filename = join(self.outputDir, 'component-%d.txt' % count)
                with open(filename, 'w') as fp2:
                    component.summarize(fp2, count)

                componentCount = len(component)
                offsets = component.offsets
                totalReads += componentCount
                print(
                    '\nConnected component %d: %d reads, covering %d offsets '
                    '(%d to %d)' % (
                        count, componentCount, len(offsets),
                        min(offsets), max(offsets)), file=fp)

                ccCounts = sorted(
                    map(len, (cc.reads
                              for cc in component.consistentComponents)),
                    reverse=True)
                if len(ccCounts) > 1:
                    print('  largest two consistent component size ratio '
                          '%.2f' % (ccCounts[0] / ccCounts[1]), file=fp)

                for j, cc in enumerate(component.consistentComponents,
                                       start=1):
                    print('  consistent sub-component %d: read count %d, '
                          'covered offset count %d.' %
                          (j, len(cc.reads), len(cc.nucleotides)), file=fp)

            print('\nIn total, %d reads were assigned to components.' %
                  totalReads, file=fp)

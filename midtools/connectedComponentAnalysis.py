from __future__ import print_function, division

from os.path import join
from collections import defaultdict, Counter

from pysam import AlignmentFile

from dark.dna import compareDNAReads
from dark.reads import Read, Reads
from dark.fasta import FastaReads
from dark.sam import samfile

from midtools.analysis import ReadAnalysis
from midtools.data import gatherData, findSignificantOffsets
from midtools.match import matchToString
from midtools.plotting import plotBaseFrequencies, plotConsistentComponents
from midtools.utils import (
    baseCountsToStr, nucleotidesToStr, commonest, fastaIdentityTable, s,
    commas, alignmentQuality)


def connectedComponentsByOffset(significantReads, threshold):
    """
    Yield sets of reads that are connected according to what significant
    offsets they cover (the nucleotides at those offsets are irrelevant at
    this point).

    @param significantReads: A C{set} of C{AlignedRead} instances, all of
        which cover at least one significant offset.
    @param threshold: A C{float} indicating what fraction of read's nucleotides
        must be identical (to those already in the component) for it to be
        allowed to join a growing component.
    @return: A generator that yields C{ComponentByOffsets} instances.
    """
    while significantReads:
        element = sorted(significantReads)[0]
        significantReads.remove(element)
        component = {element}
        offsets = set(element.significantOffsets)
        addedSomething = True
        while addedSomething:
            addedSomething = False
            reads = set()
            for read in significantReads:
                if offsets.intersection(read.significantOffsets):
                    addedSomething = True
                    reads.add(read)
                    offsets.update(read.significantOffsets)
            if reads:
                significantReads.difference_update(reads)
                component.update(reads)
        yield ComponentByOffsets(component, offsets, threshold)


class ConsistentComponent(object):
    """
    Hold information about a set of reads that share significant offsets
    and which (largely) agree on the nucleotides present at those offsets.
    """
    def __init__(self, reads, nucleotides):
        self.reads = reads
        self.nucleotides = nucleotides

    def __len__(self):
        return len(self.reads)

    def saveFasta(self, fp):
        """
        Save all reads as FASTA.

        @param fp: A file pointer to write to.
        """
        for read in sorted(self.reads):
            print(read.toString('fasta'), end='', file=fp)

    def savePaddedFasta(self, fp):
        """
        Save all reads as FASTA, padded with gaps to preserve alignment.

        @param fp: A file pointer to write to.
        """
        for read in sorted(self.reads):
            print(read.toPaddedString(), end='', file=fp)

    def consensusSequence(self, componentOffsets, infoFp):
        """
        Get a consensus sequence.

        @param componentOffsets: The C{set} of offsets in this component. This
            is *not* the same as the offsets in this consistent component
            because this consistent component may not have reads for all
            offsets.
        @param infoFp: A file pointer to write draw (and other) info to.
        @return: A C{str} consensus sequence.
        """
        sequence = []
        for offset in sorted(componentOffsets):
            if offset in self.nucleotides:
                base = commonest(
                    self.nucleotides[offset], infoFp,
                    'WARNING: consensus draw at offset %d' % offset +
                    ' %(baseCounts)s.')
            else:
                base = '-'
            sequence.append(base)
        return ''.join(sequence)

    def saveConsensus(self, count, componentOffsets, consensusFp, infoFp):
        """
        Save a consensus as FASTA.

        @param count: The C{int} number of this consistent component within
            its overall connected component.
        @param componentOffsets: The C{set} of offsets in this component. This
            is *not* the same as the offsets in this consistent component
            because this consistent component may not have reads for all
            offsets.
        @param consensusFp: A file pointer to write the consensus to.
        @param drawFp: A file pointer to write draw (and other) info to.
        """
        print(
            Read('consistent-component-%d-consensus (based on %d reads)' %
                 (count, len(self.reads)),
                 self.consensusSequence(
                     componentOffsets, infoFp)).toString('fasta'),
            file=consensusFp, end='')

    def summarize(self, fp, count, componentOffsets):
        plural = s(len(self.reads))
        print('    Consistent component %d: %d read%s, covering %d offset%s' %
              (count, len(self.reads), plural, len(self.nucleotides),
               s(len(self.nucleotides))), file=fp)
        print('    Nucleotide counts for each offset:', file=fp)
        print(nucleotidesToStr(self.nucleotides, '      '), file=fp)
        print('    Consensus sequence: %s' %
              self.consensusSequence(componentOffsets, fp), file=fp)
        print('    Read%s:' % plural, file=fp)
        for read in sorted(self.reads):
            print('     ', read, file=fp)


class ComponentByOffsets(object):
    """
    Hold information about a set of reads that share significant offsets
    regardless of the nucleotides present at those offsets. Create a list
    of subsets of these reads (ConsistentComponent instances) that are
    consistent in the nucleotides at their offsets.
    """
    def __init__(self, reads, offsets, threshold):
        self.reads = reads
        self.offsets = offsets
        self.threshold = threshold
        self.consistentComponents = list(self.findConsistentComponents())
        self._check()

    def __len__(self):
        return len(self.reads)

    def _check(self):
        selfReads = len(self)
        ccReads = sum(map(len, self.consistentComponents))
        assert selfReads == ccReads, '%d != %d' % (selfReads, ccReads)

    def summarize(self, fp, count):
        ccLengths = ', '.join(
            str(l) for l in map(len, self.consistentComponents))
        print('component %d: %d reads, covering %d offsets, split into %d '
              'consistent sub-components of lengths %s.' % (
                  count, len(self), len(self.offsets),
                  len(self.consistentComponents), ccLengths), file=fp)
        print('  offsets:', commas(self.offsets), file=fp)
        for read in sorted(self.reads):
            print('  ', read, file=fp)

        for i, cc in enumerate(self.consistentComponents, start=1):
            print(file=fp)
            cc.summarize(fp, i, self.offsets)

    def saveFasta(self, outputDir, count, verbose):
        for i, cc in enumerate(self.consistentComponents, start=1):
            filename = join(outputDir, 'component-%d-%d.fasta' % (count, i))
            if verbose > 1:
                print('      Saving component %d %d FASTA to' % (count, i),
                      filename)
            with open(filename, 'w') as fp:
                cc.saveFasta(fp)

            filename = join(outputDir, 'component-%d-%d-padded.fasta' %
                            (count, i))
            if verbose > 1:
                print('      Saving component %d %d padded FASTA to' %
                      (count, i), filename)
            with open(filename, 'w') as fp:
                cc.savePaddedFasta(fp)

    def findConsistentComponents(self):
        """
        Find sets of reads that are consistent (up to the difference threshold
        in self.threshold) according to what nucleotides they have at their
        significant offsets.
        """

        threshold = self.threshold

        # Take a copy of self.reads so we don't empty it.
        componentReads = set(self.reads)

        while componentReads:
            reads = sorted(componentReads)
            read0 = reads[0]
            ccReads = {read0}
            nucleotides = defaultdict(Counter)
            for offset in read0.significantOffsets:
                nucleotides[offset][read0.base(offset)] += 1
            rejected = set()

            # First phase:
            #
            # Add all the reads that agree exactly at all the offsets
            # covered by the set of reads so far.
            for read in reads[1:]:
                nucleotidesIfAccepted = []
                for offset in read.significantOffsets:
                    # Sanity check:
                    base = read.base(offset)
                    if offset in nucleotides:
                        # Sanity check: in phase one there can only be one
                        # nucleotide at each offset.
                        assert len(nucleotides[offset]) == 1
                        if base not in nucleotides[offset]:
                            # Not an exact match. Reject this read for now.
                            rejected.add(read)
                            break
                    nucleotidesIfAccepted.append((offset, base))
                else:
                    # We didn't break, so add this read and its nucleotides.
                    for offset, base in nucleotidesIfAccepted:
                        nucleotides[offset][base] += 1
                    reads.remove(read)
                    ccReads.add(read)

            # Second phase, part 1.
            #
            # Add in the first-round rejects that have a high enough threshold.
            # We do this before we add the bases of the rejects to nucleotides
            # because we don't want to pollute 'nucleotides' with a bunch of
            # bases from first-round rejects (because their bases could
            # overwhelm the bases of the first round acceptees, lead to the
            # acceptance of reads that would otherwise be excluded, and later
            # distort the consensus made from this component).
            acceptedRejects = set()
            offsets = set(nucleotides)
            for read in rejected:
                # Only add rejected reads with an offset overlap with the set
                # of reads already selected.
                if set(read.significantOffsets) & offsets:
                    total = matching = 0
                    for offset in read.significantOffsets:
                        if offset in nucleotides:
                            total += 1
                            if read.base(offset) in nucleotides[offset]:
                                matching += 1
                    if total and matching / total >= threshold:
                        # Looks good!
                        acceptedRejects.add(read)
                        # TODO: the next line isn't needed?
                        reads.remove(read)
                        ccReads.add(read)

            # Second phase, part 2.
            #
            # Add the nucleotides of the second-round acceptances.
            for accepted in acceptedRejects:
                for offset in accepted.significantOffsets:
                    nucleotides[offset][accepted.base(offset)] += 1

            componentReads.difference_update(ccReads)
            yield ConsistentComponent(ccReads, nucleotides)

    def saveConsensuses(self, outputDir, count, verbose):
        consensusFilename = join(
            outputDir, 'component-%d-consensuses.fasta' % count)
        infoFilename = join(
            outputDir, 'component-%d-consensuses.txt' % count)
        if verbose:
            print('      Saving component %d consensus FASTA to %s\n'
                  '      Saving component %d consensus info to %s' %
                  (count, consensusFilename, count, infoFilename))
        with open(consensusFilename, 'w') as consensusFp, open(
                infoFilename, 'w') as infoFp:
            # First write the reference sequence for this component.
            (reference,) = list(
                FastaReads(join(outputDir,
                                'reference-component-%d.fasta' % count)))
            print(reference.toString('fasta'), file=consensusFp, end='')
            for i, cc in enumerate(self.consistentComponents, start=1):
                cc.saveConsensus(i, self.offsets, consensusFp, infoFp)

        # Write out an HTML table showing the identity between the various
        # component consensuses.
        identityTableFilename = join(
            outputDir, 'component-%d-consensuses-identity.html' % count)
        if verbose:
            print('      Saving component %d consensus identity table to %s' %
                  (count, identityTableFilename))

        fastaIdentityTable(consensusFilename, identityTableFilename, verbose)


class ConnectedComponentAnalysis(ReadAnalysis):
    """
    Perform a connected component read alignment analysis for multiple
    infection detection.

    @param alignmentFiles: A C{list} of C{str} names of SAM/BAM alignment
        files.
    @param referenceGenomeFiles: A C{list} of C{str} names of FASTA files
        containing reference genomes.
    @param referenceIds: The C{str} sequence ids whose alignment should be
        analyzed. All ids must be present in the C{referenceGenomes} files.
        One of the SAM/BAM files given using C{alignmentFiles} should have an
        alignment against the given argument. If omitted, all references that
        are aligned to in the given BAM/SAM files will be analyzed.
    @param outputDir: The C{str} directory to save result files to.
    @param minReads: The C{int} minimum number of reads that must cover a
        location for it to be considered significant.
    @param homogeneousCutoff: If the most common nucleotide at a location
        occurs more than this C{float} fraction of the time (i.e., amongst all
        reads that cover the location) then the locaion will be considered
        homogeneous and therefore uninteresting.
    @param agreementThreshold: Only reads with agreeing nucleotides at
        at least this C{float} fraction of the significant sites they have in
        common will be considered connected (this is for the second phase of
        adding reads to a component.
    @param plotSAM: If C{True} save plots of where reads lie on each reference
        genome (can be slow).
    @param saveReducedFASTA: If C{True}, write out a FASTA file of the original
        input but with just the signifcant locations.
    @param verbose: The C{int}, verbosity level. Use C{0} for no output.
    """
    DEFAULT_AGREEMENT_THRESHOLD = 0.5

    def __init__(self, alignmentFiles, referenceGenomeFiles, referenceIds=None,
                 outputDir=None, minReads=ReadAnalysis.DEFAULT_MIN_READS,
                 homogeneousCutoff=ReadAnalysis.DEFAULT_HOMOGENEOUS_CUTOFF,
                 agreementThreshold=DEFAULT_AGREEMENT_THRESHOLD,
                 saveReducedFASTA=False, plotSAM=False, verbose=0):

        ReadAnalysis.__init__(
            self, alignmentFiles, referenceGenomeFiles,
            referenceIds=referenceIds,
            outputDir=outputDir, minReads=minReads,
            homogeneousCutoff=homogeneousCutoff, plotSAM=plotSAM,
            saveReducedFASTA=saveReducedFASTA, verbose=verbose)

        self.agreementThreshold = agreementThreshold

    def run(self):
        """
        Perform a read analysis for all reference sequences.
        """
        results = ReadAnalysis.run(self)

        for alignmentFile in self.alignmentFiles:
            shortAlignmentFilename = self.shortAlignmentFilename[alignmentFile]
            alignmentOutputDir = join(self.outputDir, shortAlignmentFilename)

            self._writeAlignmentHTMLSummary(results[alignmentFile],
                                            alignmentOutputDir)

        return results

    def analyzeReferenceId(self, referenceId, alignmentFile, outputDir):
        """
        Analyze the given reference id in the given alignment file (if an
        alignment to the reference id is present).

        @param referenceId: The C{str} id of the reference sequence to analyze.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the output directory.
        @return: C{None} if C{referenceId} is not present in C{alignmentFile}
            or if no significant offsets are found. Else, a C{dict} containing
            the signifcant offsets and the consensus sequence that best matches
            C{referenceId}.
        """
        analysis = self.initialReferenceIdAnalysis(
             referenceId, alignmentFile, outputDir)

        if analysis:
            (genomeLength, alignedReads, readCountAtOffset,
             baseCountAtOffset, readsAtOffset, significantOffsets,
             samFilter, paddedSAM) = analysis
        else:
            return

        components = self.findConnectedComponents(alignedReads,
                                                  significantOffsets)
        self.saveComponentFasta(components, outputDir)

        self.summarize(alignedReads, significantOffsets, components,
                       genomeLength, outputDir)

        self.saveReferenceComponents(referenceId, components, outputDir)

        self.saveComponentConsensuses(components, outputDir)

        (consensusRead, bestCcIndices, unwantedReads, wantedCcReadCount,
         consensusReadCountAtOffset,
         consensusWantedReadsBaseCountAtOffset) = (
             self.saveClosestReferenceConsensus(
                 referenceId, components, baseCountAtOffset, genomeLength,
                 alignedReads, paddedSAM.referenceInsertions, outputDir))

        unwantedCount, unalignedCount = self.saveNonConsensusReads(
            unwantedReads, alignmentFile, referenceId, outputDir)

        # Sanity check.
        if (wantedCcReadCount + unwantedCount + unalignedCount ==
                len(samFilter.queryIds)):
            self.report(
                '    All alignment file reads accounted for: '
                'wantedCcReadCount (%d) + unwantedCount (%d) + '
                'unalignedCount (%d) == SAM query count (%d)' %
                (wantedCcReadCount, unwantedCount, unalignedCount,
                 len(samFilter.queryIds)))
        else:
            raise ValueError(
                'Not all alignment file reads accounted for: '
                'wantedCcReadCount (%d) + unwantedCount (%d) + '
                'unalignedCount (%d) != SAM query count (%d)' %
                (wantedCcReadCount, unwantedCount, unalignedCount,
                 len(samFilter.queryIds)))

        self.plotConsistentComponents(
            referenceId, alignmentFile, genomeLength, components,
            significantOffsets, outputDir)

        self.saveConsensusSAM(
            alignmentFile, set(alignedReads) - unwantedReads, outputDir)

        self.saveConsensusBaseFrequencyPlot(
            referenceId, genomeLength, consensusWantedReadsBaseCountAtOffset,
            consensusReadCountAtOffset, outputDir)

        self.saveBestNonReferenceConsensus(
            referenceId, components, baseCountAtOffset, genomeLength,
            alignedReads, paddedSAM.referenceInsertions, bestCcIndices,
            outputDir)

        return {
            'consensusRead': consensusRead,
            'components': components,
            'significantOffsets': significantOffsets,
        }

    def findConnectedComponents(self, alignedReads, significantOffsets):
        """
        Find all connected components.

        @param alignedReads: A list of C{AlignedRead} instances.
        @param significantOffsets: A C{set} of signifcant offsets.
        @return: A C{list} of C{connectedComponentsByOffset} instances.
        """
        significantReads = set(read for read in alignedReads
                               if read.significantOffsets)
        components = []
        for count, component in enumerate(
                connectedComponentsByOffset(significantReads,
                                            self.agreementThreshold),
                start=1):
            components.append(component)

        # Sanity check: The significantReads set should be be empty
        # following the above processing.
        assert len(significantReads) == 0
        return components

    def saveComponentFasta(self, components, outputDir):
        """
        Save FASTA for each component.

        @param outputDir: A C{str} directory path.
        """
        self.report('    Saving component FASTA')
        for count, component in enumerate(components, start=1):
            component.saveFasta(outputDir, count, self.verbose)

    def saveReferenceComponents(self, referenceId, components, outputDir):
        """
        Save a FASTA file for the reference containing just the offsets for
        all connected components.

        @param referenceId: The C{str} id of the reference sequence.
        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param outputDir: A C{str} directory path.
        """
        reference = self.referenceGenomes[referenceId]
        for count, component in enumerate(components, start=1):
            filename = join(outputDir, 'reference-component-%d.fasta' % count)
            self.report('    Saving reference component %d to %s' %
                        (count, filename))
            read = Read(reference.id + '-component-%d' % count,
                        reference.sequence)

            Reads([read]).filter(keepSites=component.offsets).save(filename)

    def plotConsistentComponents(
            self, referenceId, alignmentFile, genomeLength, components,
            significantOffsets, outputDir):
        """
        Make a plot of all consistent connected components.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFile: The C{str} name of an alignment file.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'consistent-components-plot.html')
        self.report('    Plotting consistent connected components.')
        infoFilename = join(outputDir, 'consistent-components-plot.txt')
        plotConsistentComponents(
            referenceId, genomeLength, components, significantOffsets,
            filename, infoFilename, outputDir, titleFontSize=17,
            axisFontSize=15,
            title='%s consistent connected components<br>from file %s' %
            (referenceId, self.shortAlignmentFilename[alignmentFile]))

    def saveClosestReferenceConsensus(
            self, referenceId, components, baseCountAtOffset, genomeLength,
            alignedReads, referenceInsertions, outputDir):
        """
        Calculate and save the best consensus to a reference genome.

        @param referenceId: The C{str} id of the reference sequence.
        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param alignedReads: A list of C{AlignedRead} instances.
        @param referenceInsertions: A C{dict} keyed by read id (the read
            that would cause a reference insertion). The values are lists
            of 2-tuples, with each 2-tuple containing an offset into the
            reference sequence and the C{str} of nucleotide that would be
            inserted starting at that offset.
        @param outputDir: A C{str} directory path.
        @return: A 2-tuple with 1) the C{dark.Read} instance with the closest
            consensus to the reference, and 2) a C{list} of the best
            consistent connected components used to make the consensus.
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
            referenceSequence = reference.sequence
            nucleotides = cc.nucleotides
            count = 0
            for offset in nucleotides:
                message = drawMessage + (
                    ' offset %d: base counts' % offset) + ' %(baseCounts)s.'
                referenceBase = referenceSequence[offset]
                componentBase = commonest(nucleotides[offset], referenceBase,
                                          drawFp, message)
                count += int(componentBase == referenceBase)
            return count

        def bestConsistentComponent(component, reference, fp):
            """
            Find the consistent component in the given C{ComponentByOffsets}
            instance that best matches the passed reference sequence.

            @param component: A C{ComponentByOffsets} instance.
            @param reference: A C{Read} instance.
            @param fp: A file pointer to write information to.
            @return: The C{int} index of the best consistent component.
            """
            bestScore = -1
            bestIndex = None
            offsetCount = len(component.offsets)
            for index, cc in enumerate(component.consistentComponents):
                # To compute how good each consistent component of a
                # ComponentByOffsets instance is, it's not enough to just
                # count the matches (or the fraction of matches) in the
                # consistent component, because those components can be very
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
                    '    Consistent component %d base draw' % (index + 1))
                print('  Consistent component %d (%d reads) has %d exact '
                      'matches with the reference, out of the %d offsets it '
                      'covers (%.2f%%).'
                      % (index + 1, len(cc.reads), matchCount,
                         len(cc.nucleotides),
                         matchCount / len(cc.nucleotides) * 100.0),
                      file=fp)
                score = matchCount / offsetCount
                if score == bestScore:
                    print('    WARNING: Consistent component %d has a score '
                          '(%.2f) draw with consistent component %d' %
                          (index + 1, score, bestIndex + 1), file=fp)
                elif score > bestScore:
                    bestScore = score
                    bestIndex = index

            print('  The best consistent component is number %d.' %
                  (bestIndex + 1), file=fp)

            return bestIndex

        reference = self.referenceGenomes[referenceId]
        fields = reference.id.split(maxsplit=1)
        if len(fields) == 1:
            referenceIdRest = ''
        else:
            referenceIdRest = ' ' + fields[1]

        infoFile = join(outputDir, 'consensus.txt')
        self.report('    Saving closest consensus to reference info to',
                    infoFile)

        with open(infoFile, 'w') as infoFp:
            offsetsDone = set()
            consensus = [None] * genomeLength
            bestCcIndices = []
            for count, component in enumerate(components, start=1):
                print('\nExamining component %d with %d offsets: %s' %
                      (count, len(component.offsets),
                       commas(component.offsets)), file=infoFp)
                bestCcIndex = bestConsistentComponent(component, reference,
                                                      infoFp)
                bestCcIndices.append(bestCcIndex)
                bestCc = component.consistentComponents[bestCcIndex]
                print('  Adding best nucleotides to consensus:',
                      file=infoFp)
                for offset in sorted(bestCc.nucleotides):
                    assert consensus[offset] is None
                    referenceBase = reference.sequence[offset]
                    base = commonest(
                        bestCc.nucleotides[offset], referenceBase, infoFp,
                        ('      WARNING: base count draw at offset %d ' %
                         offset) + ' %(baseCounts)s.')
                    consensus[offset] = base
                    offsetsDone.add(offset)

                    # Do some reporting on the base just added.
                    if base == referenceBase:
                        mismatch = ''
                    else:
                        consensusBase = commonest(
                            baseCountAtOffset[offset], referenceBase, infoFp,
                            ('      WARNING: consensus base count draw at '
                             'offset %d ' % offset) + ' %(baseCounts)s.')
                        mismatch = (
                            ' (mismatch: reference has %s, all-read '
                            'consensus has %s)' % (referenceBase,
                                                   consensusBase))

                    print('    Offset %d: %s from nucleotides %s%s' %
                          (offset, base,
                           baseCountsToStr(bestCc.nucleotides[offset]),
                           mismatch), file=infoFp)

            # Make two sets of reads: 1) of all the reads in the wanted
            # consistent components, and 2) all the reads in the unwanted
            # consistent components. Do this so that we do not look at the
            # unwanted reads when filling in consensus bases from the
            # non-significant offsets.
            wantedCcReads = set()
            unwantedCcReads = set()
            for bestCcIndex, component in zip(bestCcIndices, components):
                for index, cc in enumerate(component.consistentComponents):
                    if index == bestCcIndex:
                        wantedCcReads |= cc.reads
                    else:
                        # Sanity check.
                        assert not (unwantedCcReads & cc.reads)
                        unwantedCcReads |= cc.reads

            # Get the base counts at each offset, from the full set of
            # aligned reads minus the reads we don't want because they're
            # in a consistent component that is not the best for this
            # reference sequence.
            (consensusReadCountAtOffset,
             consensusWantedReadsBaseCountAtOffset, _) = gatherData(
                 genomeLength, set(alignedReads) - unwantedCcReads)

            depthFile = join(outputDir, 'consensus-depth.txt')
            self.report('    Writing consensus depth information to',
                        depthFile)
            with open(depthFile, 'w') as depthFp:
                for offset in range(genomeLength):
                    print(offset + 1, consensusReadCountAtOffset[offset],
                          file=depthFp)

            # Fill in (from the overall read consensus) the offsets that
            # were not significant in any connected component, based only
            # on reads that were in the chosen consistent components.
            offsetsToTry = sorted(set(range(genomeLength)) - offsetsDone)
            print('\nAdding bases from %d non-connected-component '
                  'consensus offsets, EXCLUDING reads belonging to '
                  'non-optimal consistent components:' % len(offsetsToTry),
                  file=infoFp)
            for offset in offsetsToTry:
                assert consensus[offset] is None
                baseCount = consensusWantedReadsBaseCountAtOffset[offset]
                if baseCount:
                    referenceBase = reference.sequence[offset]
                    base = commonest(
                        baseCount, referenceBase, infoFp,
                        ('    WARNING: consensus base count draw at '
                         'offset %d' % offset) + ' %(baseCounts)s.')
                    print('  Offset %d: %s from nucleotides %s' %
                          (offset, base, baseCountsToStr(baseCount)),
                          file=infoFp, end='')

                    if base == referenceBase:
                        print(file=infoFp)
                    else:
                        print(' (mismatch: reference has %s)' % referenceBase,
                              file=infoFp)
                    consensus[offset] = base
                    offsetsDone.add(offset)

            # Fill in (from the overall read consensus) the offsets that
            # were not significant in any connected component, including
            # from reads that were NOT in the chosen consistent components.
            # This is the best we can do with these remaining offsets (as
            # opposed to getting gaps).
            offsetsToTry = sorted(set(range(genomeLength)) - offsetsDone)
            print('\nAdding bases from %d non-connected-component '
                  'consensus offsets, INCLUDING from reads belonging to '
                  'non-optimal consistent components:' % len(offsetsToTry),
                  file=infoFp)
            for offset in offsetsToTry:
                assert consensus[offset] is None
                referenceBase = reference.sequence[offset]
                baseCount = baseCountAtOffset[offset]
                if baseCount:
                    base = commonest(
                        baseCount, referenceBase, infoFp,
                        ('    WARNING: consensus base count draw at '
                         'offset %d' % offset) + ' %(baseCounts)s.')
                    print('  Offset %d: %s from nucleotides %s' %
                          (offset, base, baseCountsToStr(baseCount)),
                          file=infoFp, end='')
                else:
                    # The reads did not cover this offset.
                    base = '-'
                    print('  Offset %d: -' % offset, file=infoFp, end='')

                if base == referenceBase:
                    print(file=infoFp)
                else:
                    print(' (mismatch: reference has %s)' % referenceBase,
                          file=infoFp)
                consensus[offset] = base
                offsetsDone.add(offset)

            # Sanity check: make sure we processed all offsets.
            assert offsetsDone == set(range(genomeLength))

            consensusId = (
                '%s-consensus-ccc best-consistent-components:%s%s' %
                (self.shortReferenceId[referenceId],
                 ','.join(map(str, bestCcIndices)), referenceIdRest))

            consensus = Read(consensusId, ''.join(consensus))

            # Print details of the match of the consensus to the reference.
            match = compareDNAReads(reference, consensus)
            print('\nOVERALL match with reference:', file=infoFp)
            print(matchToString(match, reference, consensus, indent='  '),
                  file=infoFp)

            # Print any insertions to the reference.
            wantedReadsWithInsertions = (
                set(referenceInsertions) &
                (set(alignedReads) - unwantedCcReads))
            if wantedReadsWithInsertions:
                print('\nReference insertions present in %d read%s:' % (
                    len(wantedReadsWithInsertions),
                    s(len(wantedReadsWithInsertions))), file=infoFp)
                nucleotides = defaultdict(Counter)
                for readId in wantedReadsWithInsertions:
                    for (offset, sequence) in referenceInsertions[readId]:
                        for index, base in enumerate(sequence):
                            nucleotides[offset + index][base] += 1
                print(nucleotidesToStr(nucleotides, prefix='  '), file=infoFp)
            else:
                print('\nReference insertions: none.', file=infoFp)

        filename = join(outputDir, 'reference-consensus.fasta')
        self.report('    Saving consensus to', filename)
        Reads([consensus]).save(filename)

        wantedCcReadCount = 0
        filename = join(outputDir, 'cc-wanted.fastq')
        with open(filename, 'w') as fp:
            for wantedCcRead in wantedCcReads:
                alignment = wantedCcRead.alignment
                if not (alignment.is_secondary or alignment.is_supplementary):
                    wantedCcReadCount += 1
                    print(
                        Read(alignment.query_name,
                             alignment.query_sequence,
                             alignmentQuality(alignment)).toString('fastq'),
                        end='', file=fp)
        self.report(
            '    Saved %d read%s wanted in consistent connected components '
            'to %s' % (wantedCcReadCount, s(wantedCcReadCount), filename))

        unwantedReads = set(alignedReads) - wantedCcReads

        return (consensus, bestCcIndices, unwantedReads, wantedCcReadCount,
                consensusReadCountAtOffset,
                consensusWantedReadsBaseCountAtOffset)

    def saveNonConsensusReads(self, unwantedReads, alignmentFile, referenceId,
                              outputDir):
        """
        Save the unwanted (those not from the best consistent connected
        components used to make the consensus) reads as FASTQ.

        @param unwantedReads: A C{set} of C{AlignedRead} instances.
        @param alignmentFile: The C{str} name of an alignment file.
        @param referenceId: The C{str} id of the reference sequence.
        @param outputDir: A C{str} directory path.
        @return: A 2-C{tuple} containing the number of reads that were
            unwanted and the number from the alignment file that were not
            aligned to the reference.
        """
        seenIds = set()

        def save(alignment, fp):
            if alignment.is_secondary or alignment.is_supplementary:
                return 0
            else:
                id_ = alignment.query_name
                if id_ in seenIds:
                    raise ValueError('Already seen %s' % id_)

                seenIds.add(id_)

                print(
                    Read(id_, alignment.query_sequence,
                         alignmentQuality(alignment)).toString('fastq'),
                    file=fp, end='')
                return 1

        filename = join(outputDir, 'non-consensus-reads.fastq')
        self.report('    Saving unwanted (non-consensus) reads to',
                    filename)

        with open(filename, 'w') as fp:
            # Write out the reads that aligned to the reference but which
            # we don't want because they were in consistent connected
            # components that weren't the best for the reference.
            unwantedCount = 0
            for unwantedRead in unwantedReads:
                unwantedCount += save(unwantedRead.alignment, fp)

            # Write out reads that were in the alignment file but which
            # didn't map to the reference. They're still of interest as
            # they may map to something else.
            unalignedCount = 0
            with samfile(alignmentFile) as sam:
                for alignment in sam.fetch():
                    if alignment.reference_name != referenceId:
                        unalignedCount += save(alignment, fp)

        self.report(
            '      Wrote %d read%s from %s that mapped to %s (but were '
            'unwanted) and %d that did not' %
            (unwantedCount, s(unwantedCount),
             alignmentFile, referenceId, unalignedCount))

        return unwantedCount, unalignedCount

    def saveBestNonReferenceConsensus(
            self, referenceId, components, baseCountAtOffset, genomeLength,
            alignedReads, referenceInsertions, referenceCcIndices, outputDir):
        """
        Calculate and save the best consensus that does not include the
        consistent components that were chosen for the consensus against the
        reference. This produces the best 'other' consensus in case there was
        a double infection and one of the viruses was the reference.

        @param referenceId: The C{str} id of the reference sequence.
        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param alignedReads: A list of C{AlignedRead} instances.
        @param referenceInsertions: A C{dict} keyed by read id (the read
            that would cause a reference insertion). The values are lists
            of 2-tuples, with each 2-tuple containing an offset into the
            reference sequence and the C{str} of nucleotide that would be
            inserted starting at that offset.
        @param referenceCcIndices: A list of C{int} indices of the best
            consistent connected components against the reference. These will
            not be used in making the best non-reference consensus.
        @param outputDir: A C{str} directory path.
        @return: A C{dark.Read} instance with the best non-reference consensus.
        """

        def bestConsistentComponent(component, referenceCcIndex, fp):
            """
            Find the consistent component in the given C{ComponentByOffsets}
            instance that's best to use as a non-reference component.

            @param component: A C{ComponentByOffsets} instance.
            @param referenceCcIndex: The C{int} index of the consistent
                component that was used to make the consensus to the reference.
                That consistent component cannot be used unless there is no
                other choice.
            @param fp: A file pointer to write information to.
            @return: The C{int} index of the best consistent component.
            """
            offsetCount = len(component.offsets)

            if len(component.consistentComponents) == 1:
                assert referenceCcIndex == 0
                cc = component.consistentComponents[0]
                print('  There is only one consistent connected component! '
                      'The non-reference consensus will be the same as the '
                      'reference consensus for this set of signifcant '
                      'offsets.', file=fp)
                print('  Consistent component 1 (%d reads) has %d offsets '
                      'of the %d offsets in the connected component (%.2f%%).'
                      % (len(cc.reads), len(cc.nucleotides),
                         offsetCount,
                         len(cc.nucleotides) / offsetCount * 100.0),
                      file=fp)
                return 0

            # The bestScore tuple will hold the fraction of the connected
            # components offsets that the best consistent component covers
            # and the number of reads in the best consistent component.
            bestScore = (0.0, 0)
            bestIndex = None

            for index, cc in enumerate(component.consistentComponents):
                if index == referenceCcIndex:
                    print('  Ignoring reference consistent component %d.' %
                          (referenceCcIndex + 1), file=fp)
                    continue
                fraction = len(cc.nucleotides) / offsetCount
                print('  Consistent component %d (%d reads) has %d offsets '
                      'of the %d offsets in the connected component (%.2f%%).'
                      % (index + 1, len(cc.reads), len(cc.nucleotides),
                         offsetCount,
                         len(cc.nucleotides) / offsetCount * 100.0),
                      file=fp)
                score = (fraction, len(cc.reads))
                if score == bestScore:
                    print('    WARNING: Consistent component %d has a score '
                          '(%.2f) and read count (%d) draw with consistent '
                          'component %d' %
                          (index + 1, fraction, score[1], bestIndex + 1),
                          file=fp)
                elif score > bestScore:
                    bestScore = score
                    bestIndex = index

            print('  The best non-reference consistent component is number '
                  '%d.' % (bestIndex + 1), file=fp)

            return bestIndex

        reference = self.referenceGenomes[referenceId]
        fields = reference.id.split(maxsplit=1)
        if len(fields) == 1:
            referenceIdRest = ''
        else:
            referenceIdRest = ' ' + fields[1]

        infoFile = join(outputDir, 'non-reference-consensus.txt')
        self.report('    Saving info on best non-reference consensus to',
                    infoFile)

        with open(infoFile, 'w') as infoFp:
            offsetsDone = set()
            consensus = [None] * genomeLength
            bestCcIndices = []
            for count, (referenceCcIndex, component) in enumerate(
                    zip(referenceCcIndices, components), start=1):
                print('\nExamining component %d with %d offsets: %s' %
                      (count, len(component.offsets),
                       commas(component.offsets)), file=infoFp)
                bestCcIndex = bestConsistentComponent(
                    component, referenceCcIndex, infoFp)
                bestCcIndices.append(bestCcIndex)
                bestCc = component.consistentComponents[bestCcIndex]
                print('  Adding best nucleotides to consensus:',
                      file=infoFp)
                for offset in sorted(bestCc.nucleotides):
                    assert consensus[offset] is None
                    referenceBase = reference.sequence[offset]
                    base = commonest(
                        bestCc.nucleotides[offset], referenceBase, infoFp,
                        ('    WARNING: base count draw at offset %d ' %
                         offset) + ' %(baseCounts)s.')
                    if base == referenceBase:
                        mismatch = ''
                    else:
                        consensusBase = commonest(
                            baseCountAtOffset[offset], referenceBase, infoFp,
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

            # Make a set of all the reads in the wanted consistent
            # components, and a set of all the reads in the unwanted
            # consistent components so that we do not look at the unwanted
            # reads when filling in consensus bases from the
            # non-significant offsets.
            wantedCcReads = set()
            unwantedCcReads = set()
            for bestCcIndex, component in zip(bestCcIndices, components):
                for index, cc in enumerate(component.consistentComponents):
                    if index == bestCcIndex:
                        wantedCcReads |= cc.reads
                    else:
                        # Sanity check.
                        assert not (unwantedCcReads & cc.reads)
                        unwantedCcReads |= cc.reads

            # Get the base counts at each offset, from the full set of
            # aligned reads minus the reads we don't want because they're
            # in a consistent component that is not the best for this
            # non-reference sequence.
            consensusReadCountAtOffset, wantedReadBaseCountAtOffset, _ = (
                gatherData(genomeLength, set(alignedReads) - unwantedCcReads))

            depthFile = join(outputDir, 'non-reference-consensus-depth.txt')
            self.report('    Writing non-reference consensus depth '
                        'information to', depthFile)
            with open(depthFile, 'w') as depthFp:
                for offset in range(genomeLength):
                    print(offset + 1, consensusReadCountAtOffset[offset],
                          file=depthFp)

            # Fill in (from the overall read consensus) the offsets that
            # were not significant in any connected component, based only
            # on reads that were in the chosen consistent components.
            offsetsToTry = sorted(set(range(genomeLength)) - offsetsDone)
            print('\nAdding bases from %d non-connected-component '
                  'consensus offsets, EXCLUDING reads belonging to '
                  'non-optimal consistent components:' % len(offsetsToTry),
                  file=infoFp)
            for offset in offsetsToTry:
                assert consensus[offset] is None
                baseCount = wantedReadBaseCountAtOffset[offset]
                if baseCount:
                    referenceBase = reference.sequence[offset]
                    base = commonest(
                        baseCount, referenceBase, infoFp,
                        ('    WARNING: consensus base count draw at '
                         'offset %d' % offset) + ' %(baseCounts)s.')
                    print('  Offset %d: %s from nucleotides %s' %
                          (offset, base, baseCountsToStr(baseCount)),
                          file=infoFp, end='')

                    if base == referenceBase:
                        print(file=infoFp)
                    else:
                        print(' (mismatch: reference has %s)' % referenceBase,
                              file=infoFp)
                    consensus[offset] = base
                    offsetsDone.add(offset)

            # Fill in (from the overall read consensus) the offsets that
            # were not significant in any connected component, including
            # from reads that were NOT in the chosen consistent components.
            # This is the best we can do with these remaining offsets (as
            # opposed to getting gaps).
            offsetsToTry = sorted(set(range(genomeLength)) - offsetsDone)
            print('\nAdding bases from %d non-connected-component '
                  'consensus offsets, INCLUDING from reads belonging to '
                  'non-optimal consistent components:' % len(offsetsToTry),
                  file=infoFp)
            for offset in offsetsToTry:
                assert consensus[offset] is None
                referenceBase = reference.sequence[offset]
                baseCount = baseCountAtOffset[offset]
                if baseCount:
                    base = commonest(
                        baseCount, referenceBase, infoFp,
                        ('    WARNING: consensus base count draw at '
                         'offset %d' % offset) + ' %(baseCounts)s.')
                    print('  Offset %d: %s from nucleotides %s' %
                          (offset, base, baseCountsToStr(baseCount)),
                          file=infoFp, end='')
                else:
                    # The reads did not cover this offset.
                    base = '-'
                    print('  Offset %d: -' % offset, file=infoFp, end='')

                if base == referenceBase:
                    print(file=infoFp)
                else:
                    print(' (mismatch: reference has %s)' % referenceBase,
                          file=infoFp)
                consensus[offset] = base
                offsetsDone.add(offset)

            # Sanity check: make sure we processed all offsets.
            assert offsetsDone == set(range(genomeLength))

            consensusId = (
                '%s-non-reference-consensus-ccc '
                'best-consistent-components:%s%s' %
                (self.shortReferenceId[referenceId],
                 ','.join(map(str, bestCcIndices)), referenceIdRest))

            consensus = Read(consensusId, ''.join(consensus))

            # Print details of the match of the non-reference consensus to
            # the reference.
            match = compareDNAReads(reference, consensus)
            print('\nOVERALL match with reference:', file=infoFp)
            print(matchToString(match, reference, consensus, indent='  '),
                  file=infoFp)

            # Print any insertions to the reference.
            wantedReadsWithInsertions = (
                set(referenceInsertions) &
                (set(alignedReads) - unwantedCcReads))
            if wantedReadsWithInsertions:
                print('\nReference insertions present in %d read%s:' % (
                    len(wantedReadsWithInsertions),
                    s(len(wantedReadsWithInsertions))), file=infoFp)
                nucleotides = defaultdict(Counter)
                for readId in wantedReadsWithInsertions:
                    for (offset, sequence) in referenceInsertions[readId]:
                        for index, base in enumerate(sequence):
                            nucleotides[offset + index][base] += 1
                print(nucleotidesToStr(nucleotides, prefix='  '), file=infoFp)
            else:
                print('\nReference insertions: none.', file=infoFp)

        filename = join(outputDir, 'non-reference-consensus.fasta')
        Reads([consensus]).save(filename)

        return consensus

    def saveConsensusSAM(self, alignmentFile, wantedReads, outputDir):
        """
        Save reads from the consensus to a SAM file.

        @param alignmentFile: The C{str} name of an alignment file.
        @param wantedReads: A C{set} of wanted C{AlignedRead} instances.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'reference-consensus.sam')
        self.report('    Writing consensus SAM to', filename)
        with samfile(alignmentFile) as sam:
            alignment = AlignmentFile(filename, mode='w', template=sam)
        save = alignment.write
        for read in wantedReads:
            save(read.alignment)

    def saveConsensusBaseFrequencyPlot(
            self, referenceId, genomeLength, baseCountAtOffset,
            readCountAtOffset, outputDir):
        """
        Make a plot of the sorted base frequencies for the consensus.

        @param referenceId: The C{str} id of the reference sequence.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param readCountAtOffset: A C{list} of C{int} counts of the total
            number of reads at each genome offset (i.e., just the sum of the
            values in C{baseCountAtOffset})
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'consensus-base-frequencies.html')
        self.report('    Writing consensus base frequency plot to', filename)

        significantOffsets = list(findSignificantOffsets(
            baseCountAtOffset, readCountAtOffset, self.minReads,
            self.homogeneousCutoff))

        plotBaseFrequencies(
            significantOffsets, baseCountAtOffset, readCountAtOffset, filename,
            title='%s consensus (length %d)' % (referenceId, genomeLength),
            minReads=self.minReads, homogeneousCutoff=self.homogeneousCutoff,
            histogram=False, show=False)

    def saveComponentConsensuses(self, components, outputDir):
        """
        Write out a component consensus sequence.

        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param outputDir: A C{str} directory path.
        """
        self.report('    Saving component consensuses')
        for count, component in enumerate(components, start=1):
            component.saveConsensuses(outputDir, count, self.verbose)

    def summarize(self, alignedReads, significantOffsets, components,
                  genomeLength, outputDir):
        """
        Write out an analysis summary.

        @param alignedReads: A C{list} of C{AlignedRead} instances.
        @param significantOffsets: A C{set} of signifcant offsets.
        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'component-summary.txt')
        self.report('    Writing analysis summary to', filename)

        with open(filename, 'w') as fp:

            print('Read %d aligned reads of length %d. '
                  'Found %d significant locations.' %
                  (len(alignedReads), genomeLength,
                   len(significantOffsets)), file=fp)

            print('Reads were assigned to %d connected components:' %
                  len(components), file=fp)

            totalReads = 0
            for count, component in enumerate(components, start=1):

                filename = join(outputDir, 'component-%d.txt' % count)
                self.report('    Writing component %d summary to' % count,
                            filename)
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

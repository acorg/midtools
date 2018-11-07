from __future__ import print_function, division

from os.path import join
from collections import defaultdict, Counter

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.reads import Read, Reads

from midtools.analysis import ReadAnalysis
from midtools.clusters import ReadClusters
from midtools.data import gatherData, findSignificantOffsets
from midtools.match import matchToString
from midtools.plotting import plotBaseFrequencies, plotConsistentComponents
from midtools.utils import (
    commas, fastaIdentityTable, s, commonest, nucleotidesToStr,
    baseCountsToStr, alignmentQuality)


def connectedComponentsByOffset(significantReads, maxClusterDist, fp,
                                verbose=0):
    """
    Yield sets of reads that are connected according to what significant
    offsets they cover (the nucleotides at those offsets are irrelevant at
    this point).

    @param significantReads: A C{set} of C{AlignedRead} instances, all of
        which cover at least one significant offset.
    @param maxClusterDist: A C{float} indicating the distance beyond which
        clustering should be aborted.
    @param fp: A file-like object to write information to.
    @param verbose: The C{int}, verbosity level. Use C{0} for no output.
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
        yield ComponentByOffsets(component, offsets, maxClusterDist, fp,
                                 verbose)


class ClusterComponent(object):
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

    def consensusSequence(self, componentOffsets, referenceSequence, infoFp):
        """
        Get a consensus sequence.

        @param componentOffsets: The C{set} of offsets in this component. This
            is *not* the same as the offsets in this consistent component
            because this consistent component may not have reads for all
            offsets.
        @param referenceSequence: The C{str} reference sequence.
        @param infoFp: A file pointer to write draw (and other) info to.
        @return: A C{str} consensus sequence.
        """
        sequence = []
        for offset in sorted(componentOffsets):
            if offset in self.nucleotides:
                referenceBase = referenceSequence[offset]
                base = commonest(
                    self.nucleotides[offset], referenceBase, drawFp=infoFp,
                    drawMessage=(
                        'WARNING: consensus draw at offset %d' % offset +
                        ' %(baseCounts)s.'))
            else:
                base = '-'
            sequence.append(base)
        return ''.join(sequence)

    def saveConsensus(self, count, componentOffsets, referenceSequence,
                      consensusFp, infoFp):
        """
        Save a consensus as FASTA.

        @param count: The C{int} number of this consistent component within
            its overall connected component.
        @param componentOffsets: The C{set} of offsets in this component. This
            is *not* the same as the offsets in this consistent component
            because this consistent component may not have reads for all
            offsets.
        @param referenceSequence: The C{str} reference sequence.
        @param consensusFp: A file pointer to write the consensus to.
        @param drawFp: A file pointer to write draw (and other) info to.
        """
        sequence = self.consensusSequence(componentOffsets, referenceSequence,
                                          infoFp)
        id_ = 'consistent-component-%d-consensus (based on %d reads)' % (
               count, len(self.reads))
        print(Read(id_, sequence).toString('fasta'), file=consensusFp, end='')

    def summarize(self, fp, count, componentOffsets, referenceSequence):
        """
        Write out a summary of this consistent component.

        @param fp: The file pointer to write to.
        @param count: The C{int} number of this component.
        @param componentOffsets: The C{set} of offsets in this component.
        @param referenceSequence: The C{str} reference sequence.
        """
        plural = s(len(self.reads))
        print('    Component %d: %d read%s, covering %d offset%s' %
              (count, len(self.reads), plural, len(self.nucleotides),
               s(len(self.nucleotides))), file=fp)
        print('    Nucleotide counts for each offset:', file=fp)
        print(nucleotidesToStr(self.nucleotides, '      '), file=fp)
        print('    Consensus sequence: %s' %
              self.consensusSequence(componentOffsets, referenceSequence, fp),
              file=fp)
        print('    Read%s:' % plural, file=fp)
        for read in sorted(self.reads):
            print('     ', read, file=fp)


class ComponentByOffsets(object):
    """
    Hold information about a set of reads that share significant offsets
    regardless of the nucleotides present at those offsets. Create a list
    of subsets of these reads (ConsistentComponent instances) that are
    consistent in the nucleotides at their offsets.

    @param reads: A C{list} of C{AlignedSegment} instances.
    @param offsets: A C{set} of significant offsets covered by C{reads}.
    @param maxClusterDist: The C{float} inter-cluster distance above which
        clusters are too different to merge.
    @param fp: A file-like object to write information to.
    @param verbose: The C{int}, verbosity level. Use C{0} for no output.
    """

    def __init__(self, reads, offsets, maxClusterDist, fp, verbose=0):
        self.reads = reads
        self.offsets = offsets
        self.maxClusterDist = maxClusterDist
        self.consistentComponents = list(self.findClusters(fp, verbose))
        self._check()

    def __len__(self):
        return len(self.reads)

    def __lt__(self, other):
        """
        When components are sorted, do so by least offset.
        """
        return min(self.offsets) < min(other.offsets)

    def _check(self):
        selfReads = len(self)
        ccReads = sum(map(len, self.consistentComponents))
        assert selfReads == ccReads, '%d != %d' % (selfReads, ccReads)

    def summarize(self, fp, count, referenceSequence):
        """
        Write out a summary of this component.

        @param fp: The file pointer to write to.
        @param count: The C{int} number of this component.
        @param referenceSequence: The C{str} reference sequence.
        """
        ccLengths = ', '.join(
            str(l) for l in map(len, self.consistentComponents))
        print('component %d: %d reads, covering %d offsets, split into %d '
              'clusters of lengths %s.' % (
                  count, len(self), len(self.offsets),
                  len(self.consistentComponents), ccLengths), file=fp)
        print('  offsets:', commas(self.offsets), file=fp)
        for read in sorted(self.reads):
            print('  ', read, file=fp)

        for i, cc in enumerate(self.consistentComponents, start=1):
            print(file=fp)
            cc.summarize(fp, i, self.offsets, referenceSequence)

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

    def findClusters(self, fp, verbose):
        """
        Find clusters of reads up to the maximum cluster distance threshold in
        self.maxClusterDist according to what nucleotides they have at their
        significant offsets.

        @param fp: A file-like object to write information to.
        @param verbose: The C{int}, verbosity level. Use C{0} for no output.
        """
        readClusters = ReadClusters()

        # Using sort here makes the clustering result deterministic.
        for read in sorted(self.reads):
            readClusters.add(read)

        for readCluster in readClusters.analyze(self.maxClusterDist, fp):
            if verbose:
                print(
                    '      Found cluster with %d reads, covering %d offsets' %
                    (len(readCluster.reads), len(readCluster.nucleotides)),
                    file=fp)
            yield ClusterComponent(readCluster.reads, readCluster.nucleotides)

    def saveConsensuses(self, outputDir, count, referenceSequence, verbose):
        """
        Write out a summary of this component.

        @param outDir: The C{str} output directory.
        @param count: The C{int} number of this component.
        @param referenceSequence: The C{str} reference sequence.
        @param verbose: The C{int}, verbosity level. Use C{0} for no output.
        """
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
                cc.saveConsensus(i, self.offsets, referenceSequence,
                                 consensusFp, infoFp)

        # Write out an HTML table showing the identity between the various
        # component consensuses.
        identityTableFilename = join(
            outputDir, 'component-%d-consensuses-identity.html' % count)
        if verbose:
            print('      Saving component %d consensus identity table to %s' %
                  (count, identityTableFilename))

        fastaIdentityTable(consensusFilename, identityTableFilename, verbose)


class ClusterAnalysis(ReadAnalysis):
    """
    Perform a clustered read alignment analysis for multiple infection
    detection.

    @param alignmentFiles: A C{list} of C{str} names of SAM/BAM alignment
        files.
    @param referenceGenomeFiles: A C{list} of C{str} names of FASTA files
        containing reference genomes.
    @param referenceIds: The C{str} sequence ids whose alignment should be
        analyzed. All ids must be present in the C{referenceGenomes} files.
        One of the SAM/BAM files given using C{alignmentFiles} should have an
        alignment against the given argument. If omitted, all references that
        are aligned to in the given BAM/SAM files will be analyzed.
    @param maxClusterDist: The C{float} inter-cluster distance above which
        clusters are too different to merge.
    @param alternateNucleotideMinFreq: The C{float} frequency that an
        alternative nucleotide (i.e., not the one chosen for the consensus)
        must have in order to be selected for the alternate consensus.
    @param outputDir: The C{str} directory to save result files to.
    @param minReads: The C{int} minimum number of reads that must cover a
        location for it to be considered significant.
    @param homogeneousCutoff: If the most common nucleotide at a location
        occurs more than this C{float} fraction of the time (i.e., amongst all
        reads that cover the location) then the locaion will be considered
        homogeneous and therefore uninteresting.
    @param plotSAM: If C{True} save plots of where reads lie on each reference
        genome (can be slow).
    @param saveReducedFASTA: If C{True}, write out a FASTA file of the original
        input but with just the signifcant locations.
    @param verbose: The C{int}, verbosity level. Use C{0} for no output.
    """
    DEFAULT_MAX_CLUSTER_DIST = 0.2
    ALTERNATE_NUCLEOTIDE_MIN_FREQ_DEF = 0.15

    def __init__(
            self, alignmentFiles, referenceGenomeFiles, referenceIds=None,
            maxClusterDist=DEFAULT_MAX_CLUSTER_DIST,
            alternateNucleotideMinFreq=ALTERNATE_NUCLEOTIDE_MIN_FREQ_DEF,
            outputDir=None, minReads=ReadAnalysis.DEFAULT_MIN_READS,
            homogeneousCutoff=ReadAnalysis.DEFAULT_HOMOGENEOUS_CUTOFF,
            plotSAM=False, saveReducedFASTA=False, verbose=0):

        ReadAnalysis.__init__(
            self, alignmentFiles, referenceGenomeFiles,
            referenceIds=referenceIds,
            outputDir=outputDir, minReads=minReads,
            homogeneousCutoff=homogeneousCutoff, plotSAM=plotSAM,
            saveReducedFASTA=saveReducedFASTA, verbose=verbose)

        self.maxClusterDist = maxClusterDist
        self.alternateNucleotideMinFreq = alternateNucleotideMinFreq

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

        components = self.findConnectedComponents(
            alignedReads, significantOffsets, outputDir)

        self.saveComponentFasta(components, outputDir)

        self.summarize(referenceId, alignedReads, significantOffsets,
                       components, genomeLength, outputDir)

        self.saveReferenceComponents(referenceId, components, outputDir)

        self.saveComponentConsensuses(referenceId, components, outputDir)

        (consensusRead,
         unwantedReads,
         wantedCcReadCount,
         consensusReadCountAtOffset,
         consensusWantedReadsBaseCountAtOffset) = (
             self.saveClosestReferenceConsensus(
                 referenceId, components, baseCountAtOffset, genomeLength,
                 alignedReads, paddedSAM.referenceInsertions, outputDir))

        self.saveAlternateConsensus(
            referenceId, consensusRead, baseCountAtOffset, readCountAtOffset,
            genomeLength, outputDir)

        self.plotConsistentComponents(
            referenceId, alignmentFile, genomeLength, components,
            significantOffsets, outputDir)

        self.saveConsensusBaseFrequencyPlot(
            referenceId, genomeLength, consensusWantedReadsBaseCountAtOffset,
            consensusReadCountAtOffset, outputDir)

        return {
            'consensusRead': consensusRead,
            'components': components,
            'significantOffsets': significantOffsets,
        }

    def findConnectedComponents(self, alignedReads, significantOffsets,
                                outputDir):
        """
        Find all connected components.

        @param alignedReads: A list of C{AlignedRead} instances.
        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: The C{str} name of the output directory.
        @return: A C{list} of C{connectedComponentsByOffset} instances,
            sorted by component (the smallest offset is used for sorting
            so this gives the components from left to right along the
            reference genome.
        """
        significantReads = set(read for read in alignedReads
                               if read.significantOffsets)
        components = []
        filename = join(outputDir, 'cluster-analysis.txt')
        self.report('    Saving clustering steps to', filename)

        with open(filename, 'w') as fp:
            for count, component in enumerate(
                    connectedComponentsByOffset(
                        significantReads, self.maxClusterDist, fp,
                        self.verbose), start=1):
                components.append(component)

        # Sanity check: The significantReads set should be be empty
        # following the above processing.
        assert len(significantReads) == 0
        return sorted(components)

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
        self.report('    Plotting consistent connected components to',
                    filename)
        infoFilename = join(outputDir, 'consistent-components-plot.txt')
        self.report('    Writing consistent connected component info to',
                    infoFilename)
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
        @return: A tuple of (consensus, unwantedReads, wantedCcReadCount,
                 wantedReadsCountAtOffset, wantedReadsBaseCountAtOffset).
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
                message = (
                    drawMessage +
                    (' location %d: base counts' % (offset + 1)) +
                    ' %(baseCounts)s.')
                referenceBase = referenceSequence[offset]
                componentBase = commonest(nucleotides[offset], referenceBase,
                                          drawFp=drawFp, drawMessage=message)
                count += int(componentBase == referenceBase)
            return count

        def sortedConsistentComponent(component, reference, fp):
            """
            Sort the consistent components in the given C{ComponentByOffsets}
            instance according to how well they match the passed reference.
            The sort order is by increasing match score, so the best
            consistent component is last.

            @param component: A C{ComponentByOffsets} instance.
            @param reference: A C{Read} instance.
            @param fp: A file pointer to write information to.
            @return: The C{int} index of the best consistent component.
            """
            result = []
            for index, cc in enumerate(component.consistentComponents):
                matchCount = ccMatchCount(
                    cc, reference, fp,
                    '    Consistent component %d base draw' % (index + 1))
                score = matchCount / len(cc.nucleotides)
                print('  Consistent component %d (%d reads) has %d exact '
                      'matches with the reference, out of the %d offsets it '
                      'covers (%.2f%%).'
                      % (index + 1, len(cc.reads), matchCount,
                         len(cc.nucleotides), score * 100.0), file=fp)
                result.append((score, len(cc.nucleotides), index, cc))

            result.sort()
            return result

        reference = self.referenceGenomes[referenceId]
        fields = reference.id.split(maxsplit=1)
        if len(fields) == 1:
            referenceIdRest = ''
        else:
            referenceIdRest = ' ' + fields[1]

        infoFile = join(outputDir, 'reference-consensus.txt')
        self.report('    Saving closest consensus to reference info to',
                    infoFile)

        with open(infoFile, 'w') as infoFp:
            print('Building consensus at significant offsets.', file=infoFp)
            consensus = [None] * genomeLength
            offsetsDone = set()
            wantedReads = set()
            unwantedReads = set()
            for count, component in enumerate(components, start=1):
                print('\nExamining component %d with %d locations: %s' %
                      (count, len(component.offsets),
                       commas(
                           map(lambda offset: offset + 1, component.offsets))),
                      file=infoFp)
                componentOffsets = set(component.offsets)
                sortedCcs = sortedConsistentComponent(component, reference,
                                                      infoFp)

                while componentOffsets - offsetsDone:
                    # The following pop call will raise an IndexError if
                    # the sorted cc list is empty. But if it's empty we
                    # shouldn't be here, because the set of included
                    # offsets should at that point include everything in
                    # this component. Having the naked pop here ensures we
                    # get an exception if this assumption is incorrect.
                    # It's like having an assert to test that we found all
                    # the component's offsets following the loop.
                    score, _, ccIndex, cc = sortedCcs.pop()

                    print('  Incorporating nucleotides from consistent '
                          'component %d (%d reads, score %.2f, covering %d '
                          'locations (%d still undecided in consensus)) to '
                          'consensus.' %
                          (ccIndex + 1, len(cc.reads), score,
                           len(cc.nucleotides),
                           len(set(cc.nucleotides) - offsetsDone)),
                          file=infoFp)

                    wantedReads |= cc.reads
                    for offset in sorted(cc.nucleotides):
                        if offset in offsetsDone:
                            continue
                        nucleotides = cc.nucleotides[offset]
                        referenceBase = reference.sequence[offset]
                        base = commonest(
                            nucleotides, referenceBase, drawFp=infoFp,
                            drawMessage=(
                                '      WARNING: base count draw at '
                                'location %d '
                                % (offset + 1)) + ' %(baseCounts)s.')
                        assert consensus[offset] is None
                        consensus[offset] = base
                        offsetsDone.add(offset)

                        # Do some reporting on the base just added.
                        if base == referenceBase:
                            mismatch = ''
                        else:
                            consensusBase = commonest(
                                baseCountAtOffset[offset], referenceBase,
                                drawFp=infoFp,
                                drawMessage=(
                                    '      WARNING: consensus base count '
                                    'draw at location %d ' % (offset + 1)) +
                                ' %(baseCounts)s.')
                            mismatch = (
                                ' (mismatch: reference has %s, all-read '
                                'consensus has %s)' % (referenceBase,
                                                       consensusBase))

                        print('    Location %d: %s from nucleotides %s%s' %
                              (offset + 1, base, nucleotides.baseCountsToStr(),
                               mismatch), file=infoFp)

                # Print info about the cccs that were not needed to cover
                # all the offsets in this cc. Reverse the list so we print
                # them in decreasing match score order.
                for score, _, ccIndex, cc in reversed(sortedCcs):
                    unwantedReads |= cc.reads
                    print('  Will NOT incorporate nucleotides from consistent '
                          'component %d (%d reads, score %.2f, covering %d '
                          'locations) to consensus.' %
                          (ccIndex + 1, len(cc.reads), score,
                           len(cc.nucleotides)), file=infoFp)

            # Get the base counts at each offset, from the full set of
            # aligned reads minus the reads in cccs we're not using.
            (wantedReadsCountAtOffset,
             wantedReadsBaseCountAtOffset, _) = gatherData(
                 genomeLength, set(alignedReads) - unwantedReads)

            # Process the insignificant offsets, based on all reads EXCEPT
            # those not used in the connected components.
            offsetsToTry = sorted(set(range(genomeLength)) - offsetsDone)
            print('\nAttempting to add bases from %d non-significant '
                  'consensus locations, using all reads, EXCEPT those '
                  'belonging to unused consistent components:' %
                  len(offsetsToTry), file=infoFp)
            for offset in offsetsToTry:
                assert consensus[offset] is None
                baseCount = wantedReadsBaseCountAtOffset[offset]
                if baseCount:
                    referenceBase = reference.sequence[offset]
                    base = commonest(
                        baseCount, referenceBase, drawFp=infoFp,
                        drawMessage=(
                            '    WARNING: consensus base count draw at '
                            'location %d' % (offset + 1)) +
                        ' %(baseCounts)s.')
                    print('  Location %d: %s from nucleotides %s' %
                          (offset + 1, base, baseCountsToStr(baseCount)),
                          file=infoFp, end='')

                    if base == referenceBase:
                        print(file=infoFp)
                    else:
                        print(' (mismatch: reference has %s)' % referenceBase,
                              file=infoFp)
                    consensus[offset] = base
                    offsetsDone.add(offset)

            # Process remaining insignificant offsets, using ALL reads
            # (i.e., including those in cccs that we wanted to avoid
            # using).  At this point, this is the best we can do with these
            # final offsets (otherwise we will get gaps - which in some
            # cases may actually might be preferable because the reference
            # sequence may not be fully covered by the actual infection
            # sequence).
            offsetsToTry = sorted(set(range(genomeLength)) - offsetsDone)
            print('\nAttempting to add bases from %d non-significant '
                  'consensus locations, using all reads, INCLUDING those '
                  'belonging to unused consistent components:' %
                  len(offsetsToTry), file=infoFp)
            for offset in offsetsToTry:
                assert consensus[offset] is None
                referenceBase = reference.sequence[offset]
                baseCount = baseCountAtOffset[offset]
                if baseCount:
                    base = commonest(
                        baseCount, referenceBase, drawFp=infoFp,
                        drawMessage=(
                            '    WARNING: consensus base count draw at '
                            'location %d' % (offset + 1)) + ' %(baseCounts)s.')
                    print('  Location %d: %s from nucleotides %s' %
                          (offset + 1, base, baseCountsToStr(baseCount)),
                          file=infoFp, end='')
                else:
                    # The reads did not cover this offset.
                    base = '-'
                    print('  Location %d: -' % (offset + 1), file=infoFp,
                          end='')

                if base == referenceBase:
                    print(file=infoFp)
                else:
                    print(' (mismatch: reference has %s)' % referenceBase,
                          file=infoFp)
                consensus[offset] = base
                offsetsDone.add(offset)

            # Sanity check: make sure we processed all offsets.
            assert offsetsDone == set(range(genomeLength))

            consensusId = ('%s-consensus%s' %
                           (self.shortReferenceId[referenceId],
                            referenceIdRest))

            consensus = Read(consensusId, ''.join(consensus))

            # Print details of the match of the consensus to the reference.
            match = compareDNAReads(reference, consensus)
            print('\nOVERALL match with reference:', file=infoFp)
            print(matchToString(match, reference, consensus, indent='  '),
                  file=infoFp)

            # Print any insertions to the reference.
            wantedReadsWithInsertions = (
                set(referenceInsertions) &
                (set(alignedReads) - unwantedReads))
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
            for wantedCcRead in wantedReads:
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

        unwantedReads = set(alignedReads) - wantedReads

        return (consensus, unwantedReads, wantedCcReadCount,
                wantedReadsCountAtOffset, wantedReadsBaseCountAtOffset)

    def saveAlternateConsensus(
            self, referenceId, consensusRead, baseCountAtOffset,
            readCountAtOffset, genomeLength, outputDir):
        """
        Calculate and save an alternate consensus to a reference genome.

        @param referenceId: The C{str} id of the reference sequence.
        @consensusRead: The C{dark.reads.Read} consensus sequence to calculate
            an alternative to.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param readCountAtOffset: A C{list} of C{int} counts of the total
            number of reads at each genome offset (i.e., just the sum of the
            values in C{baseCountAtOffset})
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param outputDir: A C{str} directory path.
        @return: An alternate consensus C{dark.reads.Read} instance.
        """
        filename = join(outputDir, 'reference-alternate-consensus.txt')
        self.report('    Writing alternate consensus info to', filename)
        alternateConsensus = []

        referenceRead = self.referenceGenomes[referenceId]
        fields = referenceRead.id.split(maxsplit=1)
        if len(fields) == 1:
            referenceIdRest = ''
        else:
            referenceIdRest = ' ' + fields[1]

        with open(filename, 'w') as infoFp:
            print('The three nucleotides shown on each line are for the '
                  'reference, the consensus (made by clustering) and the '
                  'alternate nucleotide.\nThese are followed by up to two '
                  'asterisks: the first if the alternate nucleotide does not '
                  'agree with the reference, the second if it does not agree '
                  'with the consensus.', file=infoFp)
            for offset in range(genomeLength):
                referenceBase = referenceRead.sequence[offset]
                consensusBase = consensusRead.sequence[offset]
                baseCount = baseCountAtOffset[offset]
                if baseCount:
                    if len(baseCount) == 1:
                        # Only one nucleotide was found at this location.
                        # The reference doesn't necessarily agree with the
                        # consensus here, since the aligned reads may have
                        # had a different base at this site.
                        base = consensusBase
                    else:
                        # Find the nucleotide with the highest count that
                        # is not the consensus sequence nucleotide.
                        orderedCounts = baseCount.most_common()
                        alternateBase, alternateCount = [
                            x for x in orderedCounts
                            if x[0] != consensusBase][0]

                        # Check that we found a base that's not the
                        # consensus base.
                        assert alternateBase != consensusBase

                        # If the frequency of the alternate base is high
                        # enough, go with it. Else take the base from the
                        # original consensus.
                        alternateFraction = (alternateCount /
                                             readCountAtOffset[offset])
                        if (alternateCount > 1 and
                                alternateFraction >
                                self.alternateNucleotideMinFreq):
                            base = alternateBase
                        else:
                            base = consensusBase

                    agreeWithReference = referenceBase == base
                    agreeWithConsensus = consensusBase == base
                    print('Location %d: %s %s %s %s %s nucleotides %s' %
                          (offset + 1, referenceBase, consensusBase, base,
                           ' ' if agreeWithReference else '*',
                           ' ' if agreeWithConsensus else '*',
                           baseCountsToStr(baseCount)), file=infoFp)
                else:
                    # The reads did not cover this offset.
                    base = '-'
                    print('  Location %d: -' % (offset + 1), file=infoFp)

                alternateConsensus.append(base)

            alternateConsensusId = ('%s-alternate-consensus%s' %
                                    (self.shortReferenceId[referenceId],
                                     referenceIdRest))

            alternateConsensusRead = Read(alternateConsensusId,
                                          ''.join(alternateConsensus))

            # Print details of the match of the alternate consensus to the
            # reference.
            match = compareDNAReads(referenceRead, alternateConsensusRead)
            print('\nAlternate consensus match with reference:', file=infoFp)
            print(matchToString(match, referenceRead, alternateConsensusRead,
                                indent='  '), file=infoFp)

            # Print details of the match of the alternate consensus to the
            # original consensus.
            match = compareDNAReads(consensusRead, alternateConsensusRead)
            print('\nAlternate consensus match with original consensus:',
                  file=infoFp)
            print(matchToString(match, consensusRead, alternateConsensusRead,
                                indent='  '), file=infoFp)

            # Print details of the match of the original consensus to the
            # reference.
            match = compareDNAReads(referenceRead, consensusRead)
            print('\nOriginal consensus match with reference:', file=infoFp)
            print(matchToString(match, referenceRead, consensusRead,
                                indent='  '), file=infoFp)

        filename = join(outputDir, 'reference-alternate-consensus.fasta')
        self.report('    Saving alternate consensus FASTA to', filename)
        Reads([alternateConsensusRead]).save(filename)

        return alternateConsensusRead

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

    def saveComponentConsensuses(self, referenceId, components, outputDir):
        """
        Write out a component consensus sequence.

        @param referenceId: The C{str} id of the reference sequence.
        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param outputDir: A C{str} directory path.
        """
        self.report('    Saving component consensuses')
        reference = self.referenceGenomes[referenceId]
        for count, component in enumerate(components, start=1):
            component.saveConsensuses(outputDir, count,
                                      reference.sequence, self.verbose)

    def summarize(self, referenceId, alignedReads, significantOffsets,
                  components, genomeLength, outputDir):
        """
        Write out an analysis summary.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignedReads: A C{list} of C{AlignedRead} instances.
        @param significantOffsets: A C{set} of signifcant offsets.
        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'component-summary.txt')
        self.report('    Writing analysis summary to', filename)
        reference = self.referenceGenomes[referenceId]

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
                    component.summarize(fp2, count, reference.sequence)

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

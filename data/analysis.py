import sys
from tempfile import mkdtemp
from os import unlink
from os.path import exists, join, basename
from os import mkdir
from math import log10
from collections import Counter
from pathlib import Path  # This is Python 3 only.
from itertools import chain
from pprint import pprint

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.reads import Read, Reads
from dark.sam import PaddedSAM

from data.component import connectedComponentsByOffset
from data.read import AlignedRead
from data.utils import baseCountsToStr, commonest, samfile


class ReadAnalysis(object):
    """
    Perform a read alignment analysis to help with multiple infection
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
    @param saveReducedFASTA: If C{True}, write out a FASTA file of the original
        input but with just the signifcant locations.
    @param verbose: If C{True}, print verbose output.
    """
    DEFAULT_HOMOGENEOUS_CUTOFF = 0.9
    DEFAULT_MIN_READS = 5
    DEFAULT_AGREEMENT_THRESHOLD = 0.55

    def __init__(self, alignmentFiles, referenceGenomeFiles, referenceIds=None,
                 outputDir=None, minReads=DEFAULT_MIN_READS,
                 homogeneousCutoff=DEFAULT_HOMOGENEOUS_CUTOFF,
                 agreementThreshold=DEFAULT_AGREEMENT_THRESHOLD,
                 saveReducedFASTA=False, verbose=False):

        self.alignmentFiles = alignmentFiles
        self.referenceGenomeFiles = referenceGenomeFiles
        self.referenceIds = referenceIds
        self.outputDir = outputDir
        self.minReads = minReads
        self.homogeneousCutoff = homogeneousCutoff
        self.agreementThreshold = agreementThreshold
        self.saveReducedFASTA = saveReducedFASTA
        self.verbose = verbose

    def run(self):
        """
        Perform the read analysis for all reference sequences.
        """
        outputDir = self._setupOutputDir()
        self.referenceGenomes = self._readReferenceGenomes()
        self._checkReferenceIdsHaveGenomes()
        alignedReferences = self._checkReferenceIdsAreInAlignmentFiles()

        # If we weren't told which reference ids to examine the alignments of,
        # examine everything available.
        if not self.referenceIds:
            self.referenceIds = alignedReferences
            if self.verbose:
                print('Analyzing %d reference alignment%s:\n%s' %
                      (len(alignedReferences),
                       '' if len(alignedReferences) == 1 else 's',
                       '\n'.join('  %s' % id_ for id_ in alignedReferences)))

        for alignmentFile in self.alignmentFiles:
            self._analyzeAlignmentFile(alignmentFile, outputDir)

    def _analyzeAlignmentFile(self, alignmentFilename, outputDir):
        """
        Analyze all the relevant alignments in C{filename}.

        @param alignmentFilename: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the output directory.
        """
        for referenceId in self.referenceIds:
            self._analyzeReferenceId(referenceId, alignmentFilename, outputDir)

    def _analyzeReferenceId(self, referenceId, alignmentFilename, outputDir):
        """
        Analyze the given reference id in the given alignment file (if an
        alignment to the reference id is present).

        @param referenceId: The C{str} id of the reference sequence to analyze.
        @param alignmentFilename: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the output directory.

        @return: TODO..................
        """
        (thisOutputDir,
         shortAlignmentFilename,
         shortReferenceId) = self._setupAlignmentReferenceOutputDir(
                referenceId, alignmentFilename, outputDir)

        with samfile(alignmentFilename) as sam:
            tid = sam.get_tid(referenceId)
            if tid == -1:
                return
            else:
                genomeLength = sam.lengths[tid]

        alignedReads = []
        for query in PaddedSAM(alignmentFilename).queries(
                referenceName=referenceId):
            assert len(query) == genomeLength
            alignedReads.append(AlignedRead(query.id, query.sequence))

        readCountAtOffset, baseCountAtOffset, readsAtOffset = self._gatherData(
            genomeLength, alignedReads)

        significantOffsets = list(self._findSignificantOffsets(
            baseCountAtOffset, readCountAtOffset))

        if self.verbose:
            print('Read %d aligned reads of length %d. '
                  'Found %d significant locations.' %
                  (len(alignedReads), genomeLength, len(significantOffsets)))

        if not significantOffsets:
            print('No significant locations for alignment/reference combo...')
            return

        for read in alignedReads:
            read.setSignificantOffsets(significantOffsets)

        components = self._findConnectedComponents(alignedReads,
                                                   significantOffsets)

        self.saveSignificantOffsets(significantOffsets, thisOutputDir)
        self.saveComponentFasta(components, thisOutputDir)
        if self.saveReducedFASTA:
            self.saveReducedFasta(significantOffsets, thisOutputDir)
        self.summarize(alignedReads, significantOffsets, components,
                       genomeLength, thisOutputDir)
        self.saveComponentConsensuses(components, thisOutputDir)
        self.saveClosestReferenceConsensuses(
            components, baseCountAtOffset, genomeLength, thisOutputDir)

    def _findSignificantOffsets(self, baseCountAtOffset, readCountAtOffset):
        """
        Find the genome offsets that have significant base variability.

        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param readCountAtOffset: A C{list} of C{int} counts of the total
            number of reads at each genome offset (i.e., just the sum of the
            values in C{baseCountAtOffset})
        @return: A generator that yields 0-based significant offsets.
        """
        minReads = self.minReads
        homogeneousCutoff = self.homogeneousCutoff

        for offset, (readCount, counts) in enumerate(
                zip(readCountAtOffset, baseCountAtOffset)):
            if (readCount >= minReads and
                    max(counts.values()) / readCount <= homogeneousCutoff):
                yield offset

    def _gatherData(self, genomeLength, alignedReads):
        """
        Analyze the aligned reads.

        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param alignedReads: A C{list} of C{AlignedRead} instances.
        @return: A tuple of C{list}s (readCountAtOffset, baseCountAtOffset,
            readsAtOffset), each indexed from zero to the genom length.
        """
        readCountAtOffset = []
        baseCountAtOffset = []
        readsAtOffset = []

        nucleotides = set('ACGT')

        for offset in range(genomeLength):
            reads = set()
            counts = Counter()
            for read in alignedReads:
                base = read.base(offset)
                if base in nucleotides:
                    counts[base] += 1
                    reads.add(read)
            baseCountAtOffset.append(counts)
            readCountAtOffset.append(sum(counts.values()))
            readsAtOffset.append(reads)

        return readCountAtOffset, baseCountAtOffset, readsAtOffset

    def _setupOutputDir(self):
        """
        Set up the output directory and return its path.

        @return: The C{str} path of the output directory.
        """
        if self.outputDir:
            outputDir = self.outputDir
            if exists(outputDir):
                self._removePreExistingTopLevelOutputFiles()
            else:
                mkdir(outputDir)
        else:
            outputDir = mkdtemp()
            print('Writing output files to %s' % outputDir)
        return outputDir

    def _setupAlignmentReferenceOutputDir(self, referenceId,
                                          alignmentFilename, outputDir):
        """
        Set up the output directory for a given alignment file and reference.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFilename: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the top-level output directory.

        @return: A 3-tuple with the C{str} names of the output directory and
            the short alignment and reference sub-directory names.
        """
        # Make short versions of the reference id and filename for a
        # per-alignment-file per-reference-sequence output directory.
        shortReferenceId = referenceId.split()[0]
        shortAlignmentFilename = basename(alignmentFilename).rsplit(
            '.', maxsplit=1)[0]

        directory = join(outputDir, shortAlignmentFilename)
        if not exists(directory):
            mkdir(directory)

        directory = join(directory, shortReferenceId)
        if exists(directory):
            self._removePreExistingOutputFiles(directory)
        else:
            mkdir(directory)

        return directory, shortAlignmentFilename, shortReferenceId

    def _checkReferenceIdsHaveGenomes(self):
        """
        If any reference ids were given, check that we have a reference genome
        for each.

        @raise ValueError: If any reference id is not present in the
            reference genome files.
        """
        for id_ in self.referenceIds or []:
            if id_ not in self.referenceGenomes:
                raise ValueError(
                    'Reference id %r is not present in any reference genome '
                    'file.' % id_)

    def _checkReferenceIdsAreInAlignmentFiles(self):
        """
        If any reference ids were given, check that there is at least one
        alignment file that contains alignments against the reference id.

        @raise ValueError: If any reference id is not present in any alignment
            file.
        @return: A C{set} of C{str} reference ids as found in all passed
            alignment files.
        """
        # Get the names of all references in all alignment files.
        alignedReferences = set()
        for filename in self.alignmentFiles:
            with samfile(filename) as sam:
                for i in range(sam.nreferences):
                    alignedReferences.add(sam.get_reference_name(i))

        if self.referenceIds:
            missing = set(self.referenceIds) - alignedReferences
            if missing:
                raise ValueError(
                    'Alignments against the following reference id%s are not '
                    'present in any alignment file:\n%s' %
                    ('' if len(missing) == 1 else 's',
                     '\n'.join('  %s' % id_ for id_ in sorted(missing))))

        return alignedReferences

    def _readReferenceGenomes(self):
        """
        Read reference genomes from files and check that any duplicates have
        identical sequences.

        @raise ValueError: If a reference genome is found in more than one file
            and the sequences are not identical.
        @return: A C{dict} keyed by C{str} sequence id with C{dark.Read}
            values holding reference genomes.
        """
        result = {}
        seen = {}
        for filename in self.referenceGenomeFiles:
            for read in FastaReads(filename):
                id_ = read.id
                if id_ in seen:
                    if result[id_].sequence != read.sequence:
                        raise ValueError(
                            'Reference genome id %r was found in two files '
                            '(%r and %r) but with different sequences.' %
                            (id_, seen[id_], filename))
                else:
                    seen[id_] = filename
                    result[id_] = read
        if self.verbose:
            print('Read %d reference genome%s:\n%s' % (
                len(result), '' if len(result) == 1 else 's',
                '\n'.join('  %s' % id_ for id_ in result)))
        return result

    def _findConnectedComponents(self, alignedReads, significantOffsets):
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

    def _removePreExistingTopLevelOutputFiles(self):
        """
        Remove all pre-existing files from the top-level output directory.
        """
        pass

    def _removePreExistingOutputFiles(self, directory):
        """
        Remove all pre-existing files from the output directory for a
        particular reference sequence alignment.

        @param directory: The C{str} directory to examine.
        """
        # This prevents us from doing a run that results in (say) 6
        # component files and then later doing a run that results in
        # only 5 components and erroneously thinking that
        # component-6-2.fasta etc. are from the most recent run.
        paths = list(map(str, chain(
            Path(directory).glob('component-*.fasta'),
            Path(directory).glob('component-*.txt'),
            Path(directory).glob('significant-offsets.txt'))))

        if paths:
            if self.verbose:
                print('Removing %d pre-existing output file%s from %s '
                      'directory.' %
                      (len(paths), '' if len(paths) == 1 else 's', directory))
            list(map(unlink, paths))

    def saveComponentFasta(self, components, outputDir):
        """
        Save FASTA for each component.

        @param outputDir: A C{str} directory path.
        """
        if self.verbose:
            print('Saving component FASTA')
        for count, component in enumerate(components):
            component.saveFasta(outputDir, count)

    def saveSignificantOffsets(self, significantOffsets, outputDir):
        """
        Save the significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        if self.verbose:
            print('Saving significant offsets')
        with open(join(outputDir, 'significant-offsets.txt'), 'w') as fp:
            for offset in significantOffsets:
                print(offset, file=fp)

    def saveBaseFrequencies(self, outputDir):
        """
        Save the base frequencies.

        @param outputDir: A C{str} directory path.
        """
        genomeLengthWidth = int(log10(self.genomeLength)) + 1
        nucleotides = set('ACGT')

        with open(join(outputDir, 'base-frequencies.txt'), 'w') as fp:
            for offset in range(self.genomeLength):
                counts = Counter()
                for read in self.readsAtOffset[offset]:
                    base = read.base(offset)
                    if base in nucleotides:
                        counts[base] += 1
                print('Location %*d: base counts %s' %
                      (genomeLengthWidth, offset + 1, baseCountsToStr(counts)),
                      file=fp)

    def saveClosestReferenceConsensuses(self, components, baseCountAtOffset,
                                        genomeLength, outputDir):
        """
        Calculate and save the best consensus to each reference genome.

        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param outputDir: A C{str} directory path.
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

        if self.verbose:
            print('Saving closest consensuses to references')

        for referenceName in sorted(self.referenceGenomes):
            reference = self.referenceGenomes[referenceName]
            fields = reference.id.split(maxsplit=1)
            if len(fields) == 1:
                shortReferenceId = reference.id
                referenceIdRest = ''
            else:
                shortReferenceId = fields[0]
                referenceIdRest = ' ' + fields[1]

            infoFile = join(outputDir, '%s-consensus.txt' % shortReferenceId)
            with open(infoFile, 'w') as infoFp:
                offsetsDone = set()
                consensus = [None] * genomeLength
                bestCcIndices = []
                for count, component in enumerate(components, start=1):
                    print('\nExamining component %d with %d offsets: %s' %
                          (count, len(component.offsets),
                           ', '.join(map(str, sorted(component.offsets)))),
                          file=infoFp)
                    bestCcIndex = bestConsistentComponent(component, reference,
                                                          infoFp)
                    bestCcIndices.append(bestCcIndex)
                    bestCc = component.consistentComponents[bestCcIndex]
                    print('  Adding best nucleotides to consensus:',
                          file=infoFp)
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
                consensusOffsets = set(range(genomeLength)) - offsetsDone
                print('\nAdding bases from %d non-connected-component '
                      'consensus offsets:' % len(consensusOffsets),
                      file=infoFp)
                for offset in consensusOffsets:
                    assert consensus[offset] is None
                    baseCount = baseCountAtOffset[offset]
                    if baseCount:
                        base = commonest(
                            baseCount, infoFp,
                            ('    WARNING: consensus base count draw at '
                             'offset %d' % offset) + ' %(baseCounts)s.')
                        print('  Offset %d: %s from nucleotides %s' %
                              (offset, base, baseCountsToStr(baseCount)),
                              file=infoFp, end='')
                    else:
                        # The reads did not cover this offset.
                        base = '-'
                    referenceBase = reference.sequence[offset]
                    if base == referenceBase:
                        print(file=infoFp)
                    else:
                        print(' (mismatch: reference has %s)' %
                              referenceBase, file=infoFp)
                    consensus[offset] = base

                filename = join(outputDir,
                                '%s-consensus.fasta' % shortReferenceId)
                consensusId = (
                    '%s-consensus best-consistent-components:%s%s' %
                    (shortReferenceId, ','.join(map(str, bestCcIndices)),
                     referenceIdRest))

                consensusRead = Read(consensusId, ''.join(consensus))

                match = compareDNAReads(reference, consensusRead)

                identity = (
                    match['match']['identicalMatchCount'] +
                    match['match']['ambiguousMatchCount']) / genomeLength

                print('\nIdentity of reference with consensus: identity %.2f'
                      % (identity * 100.0), file=infoFp)
                print('\nMatch info:', file=infoFp)
                pprint(match, indent=4, stream=infoFp)
                Reads([reference, consensusRead]).save(filename)

    def saveReducedFasta(self, significantOffsets, outputDir):
        """
        Write out FASTA that contains reads with bases just at the
        significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        if self.verbose:
            print('Saving reduced FASTA')

        print('Saving reduced FASTA not implemented yet')
        return

        allGaps = '-' * len(significantOffsets)

        def unwanted(read):
            return (None if read.sequence == allGaps else read)

        FastaReads(self.fastaFile).filter(
            keepSites=significantOffsets).filter(
                modifier=unwanted).save(join(outputDir, 'reduced.fasta'))

    def saveComponentConsensuses(self, components, outputDir):
        """
        Write out a component consensus sequence.

        @param components: A C{list} of C{ComponentByOffsets} instances.
        @param outputDir: A C{str} directory path.
        """
        if self.verbose:
            print('Saving component consensuses')

        for count, component in enumerate(components, start=1):
            component.saveConsensuses(outputDir, count)

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
        if self.verbose:
            print('Writing analysis summary')

        with open(join(outputDir, 'component-summary.txt'), 'w') as fp:

            print('Read %d aligned reads of length %d. '
                  'Found %d significant locations.' %
                  (len(alignedReads), genomeLength,
                   len(significantOffsets)), file=fp)

            print('Reads were assigned to %d connected components:' %
                  len(components), file=fp)

            totalReads = 0
            for count, component in enumerate(components, start=1):

                filename = join(outputDir, 'component-%d.txt' % count)
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

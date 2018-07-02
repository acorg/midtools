from tempfile import mkdtemp
from os import unlink
from os.path import exists, join, basename
from os import mkdir
from math import log10
from collections import Counter
from pathlib import Path  # This is Python 3 only.
from itertools import chain
from collections import defaultdict

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.reads import Read, Reads
from dark.sam import PaddedSAM, samfile

from data.component import connectedComponentsByOffset
from data.read import AlignedRead
from data.utils import (
    baseCountsToStr, nucleotidesToStr, commonest, fastaIdentityTable)
from data.match import matchToString


class ReadAnalysis(object):
    """
    Perform a read alignment analysis for multiple infection detection.

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
    @param verbose: The C{int}, verbosity level. Use C{0} for no output.
    """
    DEFAULT_HOMOGENEOUS_CUTOFF = 0.9
    DEFAULT_MIN_READS = 5
    DEFAULT_AGREEMENT_THRESHOLD = 0.55

    def __init__(self, alignmentFiles, referenceGenomeFiles, referenceIds=None,
                 outputDir=None, minReads=DEFAULT_MIN_READS,
                 homogeneousCutoff=DEFAULT_HOMOGENEOUS_CUTOFF,
                 agreementThreshold=DEFAULT_AGREEMENT_THRESHOLD,
                 saveReducedFASTA=False, verbose=0):

        self.alignmentFiles = alignmentFiles
        self.referenceGenomeFiles = referenceGenomeFiles
        self.outputDir = outputDir
        self.minReads = minReads
        self.homogeneousCutoff = homogeneousCutoff
        self.agreementThreshold = agreementThreshold
        self.saveReducedFASTA = saveReducedFASTA
        self.verbose = verbose
        self.shortReferenceId = {}
        self.referenceGenomes = self._readReferenceGenomes()
        alignedReferences = self._checkReferenceIdsAreInAlignmentFiles(
            referenceIds)
        # If we weren't told which reference ids to examine the alignments
        # of, examine all available references.
        self.referenceIds = referenceIds or alignedReferences
        self._checkReferenceIdsHaveGenomes()

        # Make short output file names from the given reference file names.
        self.shortAlignmentFilename = dict(
            (filename, basename(filename).rsplit('.', maxsplit=1)[0])
            for filename in alignmentFiles)

    def run(self):
        """
        Perform a read analysis for all reference sequences.
        """
        outputDir = self._setupOutputDir()
        results = defaultdict(lambda: defaultdict(dict))

        for alignmentFile in self.alignmentFiles:
            if self.verbose:
                print('Analyzing alignment file', alignmentFile)
            for referenceId in self.referenceIds:
                if self.verbose:
                    print('  Looking for reference', referenceId)
                result = self._analyzeReferenceId(referenceId, alignmentFile,
                                                  outputDir)
                if result:
                    results[alignmentFile][referenceId] = result

            self._writeAlignmentHTMLSummary(alignmentFile,
                                            results[alignmentFile], outputDir)
        self._writeOverallResultSummary(results, outputDir)

    def _writeAlignmentHTMLSummary(self, alignmentFile, result, outputDir):
        """
        Write an HTML summary of the overall results.

        @param alignmentFile: A C{str} alignment file name.
        @param result: A C{dict} keyed by C{str} short reference name, and
           with values being C{dict}s with signifcant offsets and best
           consensus sequence for the corresponding reference in the alignment
           file.
        """
        outputDir = join(outputDir, self.shortAlignmentFilename[alignmentFile])

        referencesFilename = join(outputDir, 'references.fasta')
        if self.verbose:
            print('  Writing reference FASTA to', referencesFilename)
        with open(referencesFilename, 'w') as fp:
            for referenceId in result:
                print(self.referenceGenomes[referenceId].toString('fasta'),
                      file=fp, end='')

        consensusesFilename = join(outputDir, 'consensuses.fasta')
        if self.verbose:
            print('  Writing consensus FASTA for all references to',
                  consensusesFilename)
        with open(consensusesFilename, 'w') as fp:
            for referenceId in result:
                print(result[referenceId]['consensusRead'].toString('fasta'),
                      file=fp, end='')

        htmlFilename = join(outputDir, 'consensus-vs-reference.html')
        if self.verbose:
            print('  Writing consensus vs reference identity table to',
                  htmlFilename)
        fastaIdentityTable(consensusesFilename, htmlFilename, self.verbose,
                           filename2=referencesFilename)

        htmlFilename = join(outputDir, 'consensus-vs-consensus.html')
        if self.verbose:
            print('  Writing consensus vs consensus identity table to',
                  htmlFilename)
        fastaIdentityTable(consensusesFilename, htmlFilename, self.verbose)

    def _writeOverallResultSummary(self, results, outputDir):
        """
        Write a summary of the overall results.

        @param results: A C{dict} of C{dicts}. Keyed by C{str} short alignment
           file name, then C{str} short reference name, and with values being
           C{dict}s with signifcant offsets and best consensus sequence for
           the corresponding reference in the alignment file.
        """
        filename = join(outputDir, 'result-summary.txt')
        if self.verbose:
            print('Writing overall result summary to', filename)
        with open(filename, 'w') as fp:
            for alignmentFilename in results:
                print('Alignment file', alignmentFilename, file=fp)
                for referenceId in results[alignmentFilename]:
                    result = results[alignmentFilename][referenceId]
                    referenceRead = self.referenceGenomes[referenceId]
                    consensusRead = result['consensusRead']
                    genomeLength = len(referenceRead)
                    significantOffsets = result['significantOffsets']
                    print('\n  Reference %s (length %d)' %
                          (referenceId, genomeLength), file=fp)
                    print('    %d significant offsets found.' %
                          len(significantOffsets), file=fp)
                    print('    %d connected components.' %
                          len(result['components']), file=fp)

                    # Overall match.
                    match = compareDNAReads(referenceRead, consensusRead)
                    print('\n    Overall match of reference with consensus:',
                          file=fp)
                    print(matchToString(
                        match, referenceRead, consensusRead, indent='    '),
                          file=fp)

                    # Significant sites match.
                    match = compareDNAReads(referenceRead, consensusRead,
                                            offsets=significantOffsets)
                    print('\n    Match of reference with consensus at '
                          '%d SIGNIFICANT sites:' % len(significantOffsets),
                          file=fp)
                    print(matchToString(
                        match, referenceRead, consensusRead, indent='    ',
                        offsets=significantOffsets), file=fp)

                    # Non-significant sites match.
                    nonSignificantOffsets = (set(range(genomeLength)) -
                                             set(significantOffsets))
                    match = compareDNAReads(referenceRead, consensusRead,
                                            offsets=nonSignificantOffsets)
                    print('\n    Match of reference with consensus at '
                          '%d NON-SIGNIFICANT sites:' %
                          len(nonSignificantOffsets), file=fp)
                    print(matchToString(
                        match, referenceRead, consensusRead, indent='    ',
                        offsets=nonSignificantOffsets), file=fp)

    def _analyzeReferenceId(self, referenceId, alignmentFile, outputDir):
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
        thisOutputDir = self._setupAlignmentReferenceOutputDir(
            referenceId, alignmentFile, outputDir)

        with samfile(alignmentFile) as sam:
            tid = sam.get_tid(referenceId)
            if tid == -1:
                # This referenceId is not in this alignment file.
                if self.verbose:
                    print('    Reference %s not in alignment file.' %
                          referenceId)
                return
            else:
                genomeLength = sam.lengths[tid]
                # Sanity check.
                assert genomeLength == len(self.referenceGenomes[referenceId])

        alignedReads = []
        paddedSAM = PaddedSAM(alignmentFile)
        for query in paddedSAM.queries(referenceName=referenceId):
            assert len(query) == genomeLength
            alignedReads.append(AlignedRead(query.id, query.sequence))

        readCountAtOffset, baseCountAtOffset, readsAtOffset = self._gatherData(
            genomeLength, alignedReads)

        self.saveBaseFrequencies(thisOutputDir, genomeLength,
                                 baseCountAtOffset)

        significantOffsets = list(self._findSignificantOffsets(
            baseCountAtOffset, readCountAtOffset))

        if self.verbose:
            print('    Read %d aligned reads of length %d. '
                  'Found %d significant locations.' %
                  (len(alignedReads), genomeLength, len(significantOffsets)))

        if not significantOffsets:
            if self.verbose:
                print('    No significant locations found.')
            return

        self.saveSignificantOffsets(significantOffsets, thisOutputDir)

        for read in alignedReads:
            read.setSignificantOffsets(significantOffsets)

        components = self._findConnectedComponents(alignedReads,
                                                   significantOffsets)
        self.saveComponentFasta(components, thisOutputDir)

        if self.saveReducedFASTA:
            self.saveReducedFasta(significantOffsets, thisOutputDir)

        self.summarize(alignedReads, significantOffsets, components,
                       genomeLength, thisOutputDir)

        self.saveComponentConsensuses(components, thisOutputDir)

        consensusRead, bestCcIndices = self.saveClosestReferenceConsensus(
            referenceId, components, baseCountAtOffset, genomeLength,
            alignedReads, paddedSAM.referenceInsertions, thisOutputDir)

        nonReferenceConsensusRead = self.saveBestNonReferenceConsensus(
            referenceId, components, baseCountAtOffset, genomeLength,
            alignedReads, paddedSAM.referenceInsertions, bestCcIndices,
            thisOutputDir)

        return {
            'consensusRead': consensusRead,
            'components': components,
            'nonReferenceConsensusRead': nonReferenceConsensusRead,
            'significantOffsets': significantOffsets,
        }

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
                                          alignmentFile, outputDir):
        """
        Set up the output directory for a given alignment file and reference.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the top-level output directory.

        @return: A 3-tuple with the C{str} names of the output directory and
            the short alignment and reference sub-directory names.
        """
        # Make short versions of the reference id and filename for a
        # per-alignment-file per-reference-sequence output directory.

        shortAlignmentFilename = self.shortAlignmentFilename[alignmentFile]
        shortReferenceId = self.shortReferenceId[referenceId]

        directory = join(outputDir, shortAlignmentFilename)
        if not exists(directory):
            mkdir(directory)

        directory = join(directory, shortReferenceId)
        if exists(directory):
            self._removePreExistingOutputFiles(directory)
        else:
            mkdir(directory)

        return directory

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

    def _checkReferenceIdsAreInAlignmentFiles(self, referenceIds=None):
        """
        If any reference ids were given, check that there is at least one
        alignment file that contains alignments against the reference id.

        @param referenceIds: An iterable of C{str} sequence reference ids.
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

        if referenceIds:
            missing = set(referenceIds) - alignedReferences
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
        identical sequences. Populate C{self.shortReferenceId}.

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
                    self.shortReferenceId[id_] = id_.split()[0]

        if self.verbose > 1:
            print('Read %d reference genome%s:\n%s' % (
                len(result), '' if len(result) == 1 else 's',
                '\n'.join('  %s' % id_ for id_ in sorted(result))))
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
        paths = list(map(str, chain(
            Path(self.outputDir).glob('result-summary.fasta'))))

        if paths:
            if self.verbose > 1:
                print('    Removing %d pre-existing output file%s from '
                      'top-level output directory %s.' %
                      (len(paths), '' if len(paths) == 1 else 's',
                       self.outputDir))
            list(map(unlink, paths))

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
            Path(directory).glob('*.fasta'),
            Path(directory).glob('*.html'),
            Path(directory).glob('*.txt'))))

        if paths:
            if self.verbose > 1:
                print('    Removing %d pre-existing output file%s from %s '
                      'directory.' %
                      (len(paths), '' if len(paths) == 1 else 's', directory))
            list(map(unlink, paths))

    def saveComponentFasta(self, components, outputDir):
        """
        Save FASTA for each component.

        @param outputDir: A C{str} directory path.
        """
        if self.verbose:
            print('    Saving component FASTA:')
        for count, component in enumerate(components, start=1):
            component.saveFasta(outputDir, count, self.verbose)

    def saveSignificantOffsets(self, significantOffsets, outputDir):
        """
        Save the significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'significant-offsets.txt')
        if self.verbose:
            print('    Saving significant offsets to', filename)
        with open(filename, 'w') as fp:
            for offset in significantOffsets:
                print(offset, file=fp)

    def saveBaseFrequencies(self, outputDir, genomeLength, baseCountAtOffset):
        """
        Save the base nucleotide frequencies.

        @param outputDir: A C{str} directory path.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        """
        filename = join(outputDir, 'base-frequencies.txt')
        if self.verbose:
            print('    Saving base nucleotide frequencies to', filename)

        genomeLengthWidth = int(log10(genomeLength)) + 1

        with open(filename, 'w') as fp:
            for offset in range(genomeLength):
                print(
                    'Location %*d: base counts %s' %
                    (genomeLengthWidth, offset + 1, baseCountAtOffset[offset]),
                    file=fp)

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
        if self.verbose:
            print('    Saving closest consensus to reference to', infoFile)

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
            # reference sequence.
            consensusReadCountAtOffset, wantedReadBaseCountAtOffset, _ = (
                self._gatherData(genomeLength,
                                 set(alignedReads) - unwantedCcReads))

            depthFile = join(outputDir, 'consensus-depth.txt')
            if self.verbose:
                print('    Writing consensus depth information to', depthFile)
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
                '%s-consensus best-consistent-components:%s%s' %
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
                    '' if len(wantedReadsWithInsertions) == 1 else 's'),
                      file=infoFp)
                nucleotides = defaultdict(Counter)
                for readId in wantedReadsWithInsertions:
                    for (offset, sequence) in referenceInsertions[readId]:
                        for index, base in enumerate(sequence):
                            nucleotides[offset + index][base] += 1
                print(nucleotidesToStr(nucleotides, prefix='  '), file=infoFp)
            else:
                print('\nReference insertions: none.', file=infoFp)

        filename = join(outputDir, 'consensus.fasta')
        Reads([reference, consensus]).save(filename)

        return consensus, bestCcIndices

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
        if self.verbose:
            print('    Saving best non-reference consensus to', infoFile)

        with open(infoFile, 'w') as infoFp:
            offsetsDone = set()
            consensus = [None] * genomeLength
            bestCcIndices = []
            for count, (referenceCcIndex, component) in enumerate(
                    zip(referenceCcIndices, components), start=1):
                print('\nExamining component %d with %d offsets: %s' %
                      (count, len(component.offsets),
                       ', '.join(map(str, sorted(component.offsets)))),
                      file=infoFp)
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
                self._gatherData(genomeLength,
                                 set(alignedReads) - unwantedCcReads))

            depthFile = join(outputDir, 'non-reference-consensus-depth.txt')
            if self.verbose:
                print('    Writing non-reference consensus depth information '
                      'to', depthFile)
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
                '%s-non-reference-consensus best-consistent-components:%s%s' %
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
                    '' if len(wantedReadsWithInsertions) == 1 else 's'),
                      file=infoFp)
                nucleotides = defaultdict(Counter)
                for readId in wantedReadsWithInsertions:
                    for (offset, sequence) in referenceInsertions[readId]:
                        for index, base in enumerate(sequence):
                            nucleotides[offset + index][base] += 1
                print(nucleotidesToStr(nucleotides, prefix='  '), file=infoFp)
            else:
                print('\nReference insertions: none.', file=infoFp)

        filename = join(outputDir, 'non-reference-consensus.fasta')
        Reads([reference, consensus]).save(filename)

        return consensus

    def saveReducedFasta(self, significantOffsets, outputDir):
        """
        Write out FASTA that contains reads with bases just at the
        significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        if self.verbose:
            print('    Saving reduced FASTA')

        print('    Saving reduced FASTA not implemented yet')
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
            print('    Saving component consensuses:')

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
        if self.verbose:
            print('    Writing analysis summary to', filename)

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
                if self.verbose:
                    print('    Writing component %d summary to' % count,
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

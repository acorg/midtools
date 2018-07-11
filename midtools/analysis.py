from __future__ import print_function, division

import sys
from tempfile import mkdtemp
from os import unlink
from os.path import exists, join, basename
from os import mkdir
from math import log10
from collections import Counter
from pathlib import Path  # This is Python 3 only.
from itertools import chain
from collections import defaultdict

from pysam import AlignmentFile

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Read, Reads
from dark.sam import PaddedSAM, samfile

from midtools.component import connectedComponentsByOffset
from midtools.data import gatherData, findSignificantOffsets
from midtools.plotting import (
    plotBaseFrequencies, plotCoverageAndSignificantLocations)
from midtools.read import AlignedRead
from midtools.utils import (
    baseCountsToStr, nucleotidesToStr, commonest, fastaIdentityTable, s,
    alignmentQuality, commas)
from midtools.match import matchToString


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
    DEFAULT_AGREEMENT_THRESHOLD = 0.5

    def __init__(self, alignmentFiles, referenceGenomeFiles, referenceIds=None,
                 outputDir=None, minReads=DEFAULT_MIN_READS,
                 homogeneousCutoff=DEFAULT_HOMOGENEOUS_CUTOFF,
                 agreementThreshold=DEFAULT_AGREEMENT_THRESHOLD,
                 saveReducedFASTA=False, verbose=0):

        self.alignmentFiles = alignmentFiles
        self.outputDir = outputDir
        self.minReads = minReads
        self.homogeneousCutoff = homogeneousCutoff
        self.agreementThreshold = agreementThreshold
        self.saveReducedFASTA = saveReducedFASTA
        self.verbose = verbose
        self.referenceGenomes = self._readReferenceGenomes(
            referenceGenomeFiles)

        # Make short reference ids from the reference genomes.
        self.shortReferenceId = dict(
            (id_, id_.split()[0]) for id_ in self.referenceGenomes)

        # Make short output file names from the given reference file names.
        self.shortAlignmentFilename = dict(
            (filename, basename(filename).rsplit('.', maxsplit=1)[0])
            for filename in alignmentFiles)

        alignedReferences = self._getAlignedReferences(alignmentFiles)
        self.referenceIds = self._getReferenceIds(alignedReferences,
                                                  referenceIds)

    def _getReferenceIds(self, alignedReferences, referenceIds):
        """
        Figure out which reference ids we can process.

        @param alignedReferences: A C{set} of C{str} reference ids found in
            the passed reference files.
        @param referenceIds: A C{list} of C{str} reference ids for which
            processing has specifically been requested, or C{None}.
        @return: A C{set} of C{str} reference ids to process.
        """
        if referenceIds:
            # Specific reference ids were given. Check that each appears in
            # some alignment file and that we have a genome for each. Any
            # error here causes a message to stderr and exit.
            missing = set(referenceIds) - alignedReferences
            if missing:
                print(
                    'Alignments against the following reference id%s are not '
                    'present in any alignment file:\n%s' %
                    (s(len(missing)), '\n'.join('  %s' % id_
                                                for id_ in sorted(missing))),
                    file=sys.stderr)
                sys.exit(1)

            missing = set(referenceIds) - set(self.referenceGenomes)
            if missing:
                print(
                    'Reference id%s %s not present in any reference genome '
                    'file.' % (s(len(missing)), commas(missing)),
                    file=sys.stderr)
                sys.exit(1)
        else:
            # We weren't told which reference ids to specifically examine
            # the alignments of, so examine all available references
            # mentioned in any alignment file and that we also have a
            # genome for. Mention any references from alignment files that
            # we can't process due to lack of genome.
            missing = alignedReferences - set(self.referenceGenomes)
            if missing:
                self.report(
                    'No analysis will be performed on reference%s %s '
                    '(found in SAM/BAM alignment file(s) headers) because no '
                    'corresponding reference genome was found.' %
                    (s(len(missing)), commas(missing)))

            referenceIds = alignedReferences & set(self.referenceGenomes)

            if referenceIds:
                self.report(
                    'Examining %d reference%s: %s' %
                    (len(referenceIds), s(len(referenceIds)),
                     commas(referenceIds)))
            else:
                print(
                    'Nothing to do! No genome could be found for any aligned '
                    'reference. Found reference%s: %s' %
                    (s(len(alignedReferences)), commas(alignedReferences)),
                    file=sys.stderr)
                sys.exit(1)

        return referenceIds

    def report(self, *args, requiredVerbosityLevel=1):
        """
        Print a status message, if our verbose setting is high enough.

        @param args: The arguments to print.
        @param requiredVerbosityLevel: The minimum C{int} verbosity
            level required.
        """
        if self.verbose >= requiredVerbosityLevel:
            print(*args)

    def run(self):
        """
        Perform a read analysis for all reference sequences.
        """
        outputDir = self._setupOutputDir()
        results = defaultdict(lambda: defaultdict(dict))

        for alignmentFile in self.alignmentFiles:
            self.report('Analyzing alignment file', alignmentFile)
            alignmentOutputDir = self._setupAlignmentOutputDir(alignmentFile,
                                                               outputDir)

            self._writeAlignmentFileSummary(alignmentFile, alignmentOutputDir)

            for referenceId in self.referenceIds:
                self.report('  Looking for reference', referenceId)

                referenceOutputDir = self._setupReferenceOutputDir(
                    referenceId, alignmentOutputDir)

                result = self._analyzeReferenceId(
                    referenceId, alignmentFile, referenceOutputDir)

                if result:
                    results[alignmentFile][referenceId] = result

            self._writeAlignmentHTMLSummary(results[alignmentFile],
                                            alignmentOutputDir)

        self._writeOverallResultSummary(results, outputDir)
        self._writeOverallResultSummarySummary(results, outputDir)

    def _writeAlignmentFileSummary(self, alignmentFile, outputDir):
        """
        Write a summary of alignments.

        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the output directory.
        """
        shortAlignmentFilename = self.shortAlignmentFilename[alignmentFile]
        filename = join(outputDir, shortAlignmentFilename + '.stats')
        self.report('  Writing alignment statistics to', filename)
        e = Executor()
        e.execute('sam-reference-read-counts.py "%s" > %s' %
                  (alignmentFile, filename))
        if self.verbose > 1:
            for line in e.log:
                print('    ', line)

    def _writeAlignmentHTMLSummary(self, result, outputDir):
        """
        Write an HTML summary of the overall results.

        @param result: A C{dict} keyed by C{str} short reference name, and
           with values being C{dict}s with signifcant offsets and best
           consensus sequence for the corresponding reference in the alignment
           file.
        """
        referencesFilename = join(outputDir, 'references.fasta')
        self.report('  Writing FASTA for mapped-to references to',
                    referencesFilename)
        with open(referencesFilename, 'w') as fp:
            for referenceId in sorted(result):
                print(self.referenceGenomes[referenceId].toString('fasta'),
                      file=fp, end='')

        consensusesFilename = join(outputDir, 'consensuses.fasta')
        self.report('  Writing FASTA consensus for mapped-to references to',
                    consensusesFilename)
        with open(consensusesFilename, 'w') as fp:
            for referenceId in sorted(result):
                print(result[referenceId]['consensusRead'].toString('fasta'),
                      file=fp, end='')

        htmlFilename = join(outputDir, 'consensus-vs-reference.html')
        self.report('  Writing consensus vs reference identity table to',
                    htmlFilename)
        fastaIdentityTable(consensusesFilename, htmlFilename, self.verbose,
                           filename2=referencesFilename)

        htmlFilename = join(outputDir, 'consensus-vs-consensus.html')
        self.report('  Writing consensus vs consensus identity table to',
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
        self.report('Writing overall result summary to', filename)
        with open(filename, 'w') as fp:
            for alignmentFilename in sorted(results):
                print('Alignment file', alignmentFilename, file=fp)
                for referenceId in sorted(results[alignmentFilename]):
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

    def _writeOverallResultSummarySummary(self, results, outputDir):
        """
        Write a summary of the summary of the overall results.

        @param results: A C{dict} of C{dicts}. Keyed by C{str} short alignment
           file name, then C{str} short reference name, and with values being
           C{dict}s with signifcant offsets and best consensus sequence for
           the corresponding reference in the alignment file.
        """
        filename = join(outputDir, 'result-summary-summary.txt')
        self.report('Writing overall result summary summary to', filename)

        bestFraction = 0.0
        bestAlignmentReference = []

        with open(filename, 'w') as fp:
            for alignmentFilename in sorted(results):
                print(alignmentFilename, file=fp)
                resultSummary = []
                for referenceId in sorted(results[alignmentFilename]):
                    result = results[alignmentFilename][referenceId]
                    referenceRead = self.referenceGenomes[referenceId]
                    consensusRead = result['consensusRead']
                    match = compareDNAReads(
                        referenceRead, consensusRead)['match']
                    matchCount = (match['identicalMatchCount'] +
                                  match['ambiguousMatchCount'])
                    fraction = matchCount / len(referenceRead)

                    if fraction > bestFraction:
                        bestFraction = fraction
                        bestAlignmentReference = [
                            (alignmentFilename, referenceId)]
                    elif fraction == bestFraction:
                        bestAlignmentReference.append(
                            (alignmentFilename, referenceId))

                    resultSummary.append(
                        (fraction,
                         '  %s: %d/%d (%.2f%%)' % (
                             referenceId, matchCount, len(referenceRead),
                             fraction * 100.0)))

                # Sort the result summary by decreasing nucleotide identity
                # fraction.
                resultSummary.sort(reverse=True)
                for fraction, summary in resultSummary:
                    print(summary, file=fp)

                print(file=fp)

            print('Best match%s (%.2f%%):' %
                  ('' if len(bestAlignmentReference) == 1 else 'es',
                   bestFraction * 100.0), file=fp)
            for alignmentFilename, referenceId in bestAlignmentReference:
                print('  %s: %s' % (alignmentFilename, referenceId), file=fp)

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

        # Make sure this reference id is in this alignment file and if so
        # get its length (and check it's the same as the length of the
        # sequence given in the reference file).
        with samfile(alignmentFile) as sam:
            tid = sam.get_tid(referenceId)
            if tid == -1:
                # This referenceId is not in this alignment file.
                self.report('    Reference %s not in alignment file.' %
                            referenceId)
                return
            else:
                genomeLength = sam.lengths[tid]
                # Sanity check.
                assert genomeLength == len(self.referenceGenomes[referenceId])

        alignedReads = []
        paddedSAM = PaddedSAM(alignmentFile)
        for query in paddedSAM.queries(referenceName=referenceId,
                                       addAlignment=True, storeQueryIds=True):
            assert len(query) == genomeLength
            alignedReads.append(
                AlignedRead(query.id, query.sequence, query.alignment))

        # Sanity check that all aligned reads have different ids. This
        # should be the case because the padded SAM queries method adds /2,
        # /3 etc to queries that have more than one alignment.
        assert len(alignedReads) == len(set(read.id for read in alignedReads))

        readCountAtOffset, baseCountAtOffset, readsAtOffset = gatherData(
            genomeLength, alignedReads)

        significantOffsets = list(findSignificantOffsets(
            baseCountAtOffset, readCountAtOffset, self.minReads,
            self.homogeneousCutoff))

        self.report('    %d alignment%s (of %d unique %s) read from %s' %
                    (paddedSAM.alignmentCount,
                     s(paddedSAM.alignmentCount),
                     len(paddedSAM.queryIds),
                     'query' if len(paddedSAM.queryIds) == 1 else 'queries',
                     alignmentFile))
        self.report('    %d of which %s aligned to %s' %
                    (len(alignedReads),
                     'was' if len(alignedReads) == 1 else 'were', referenceId))
        self.report('    Reference genome length %d' % genomeLength)
        self.report('    Found %d significant location%s' %
                    (len(significantOffsets), s(len(significantOffsets))))

        self.saveBaseFrequencies(outputDir, genomeLength, baseCountAtOffset)

        if not significantOffsets:
            self.report('    No significant locations found.')
            return

        self._plotCoverageAndSignificantLocations(
            referenceId, alignmentFile, readCountAtOffset, genomeLength,
            significantOffsets, outputDir)

        self.saveSignificantOffsets(significantOffsets, outputDir)

        for read in alignedReads:
            read.setSignificantOffsets(significantOffsets)

        self.saveReferenceBaseFrequencyPlot(
            referenceId, genomeLength, significantOffsets,
            baseCountAtOffset, readCountAtOffset, outputDir)

        components = self._findConnectedComponents(alignedReads,
                                                   significantOffsets)
        self.saveComponentFasta(components, outputDir)

        if self.saveReducedFASTA:
            self.saveReducedFasta(significantOffsets, outputDir)

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
                len(paddedSAM.queryIds)):
            self.report(
                '    All alignment file reads accounted for: '
                'wantedCcReadCount (%d) + unwantedCount (%d) + '
                'unalignedCount (%d) == SAM query count (%d)' %
                (wantedCcReadCount, unwantedCount, unalignedCount,
                 len(paddedSAM.queryIds)))
        else:
            raise ValueError(
                'Not all alignment file reads accounted for: '
                'wantedCcReadCount (%d) + unwantedCount (%d) + '
                'unalignedCount (%d) != SAM query count (%d)' %
                (wantedCcReadCount, unwantedCount, unalignedCount,
                 len(paddedSAM.queryIds)))

        self.saveConsensusSAM(
            alignmentFile, set(alignedReads) - unwantedReads, outputDir)

        self.saveConsensusBaseFrequencyPlot(
            referenceId, genomeLength, consensusWantedReadsBaseCountAtOffset,
            consensusReadCountAtOffset, outputDir)

        self.saveBestNonReferenceConsensus(
            referenceId, components, baseCountAtOffset, genomeLength,
            alignedReads, paddedSAM.referenceInsertions, bestCcIndices,
            outputDir)

        # Extract a consensus according to bcftools.
        self.writeBcftoolsConsensus(referenceId, alignmentFile, outputDir)

        return {
            'consensusRead': consensusRead,
            'components': components,
            'significantOffsets': significantOffsets,
        }

    def _setupOutputDir(self):
        """
        Set up the output directory and return its path.

        @return: The C{str} path of the output directory.
        """
        if self.outputDir:
            outputDir = self.outputDir
            if exists(outputDir):
                self._removePreExistingTopLevelOutputDirFiles()
            else:
                mkdir(outputDir)
        else:
            outputDir = mkdtemp()
            print('Writing output files to %s' % outputDir)
        return outputDir

    def _setupAlignmentOutputDir(self, alignmentFile, outputDir):
        """
        Set up the output directory for a given alignment file.

        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the top-level output directory.
        @return: The C{str} output directory name.
        """
        shortAlignmentFilename = self.shortAlignmentFilename[alignmentFile]

        directory = join(outputDir, shortAlignmentFilename)
        if exists(directory):
            self._removePreExistingAlignmentDirFiles(directory)
        else:
            mkdir(directory)

        return directory

    def _setupReferenceOutputDir(self, referenceId, outputDir):
        """
        Set up the output directory for a given alignment file and reference.

        @param referenceId: The C{str} id of the reference sequence.
        @param outputDir: The C{str} name of the top-level output directory.
        @return: The C{str} output directory name.
        """
        # Make short versions of the reference id and filename for a
        # per-alignment-file per-reference-sequence output directory.

        shortReferenceId = self.shortReferenceId[referenceId]
        directory = join(outputDir, shortReferenceId)
        if exists(directory):
            self._removePreExistingReferenceDirFiles(directory)
        else:
            mkdir(directory)

        return directory

    def _getAlignedReferences(self, alignmentFiles):
        """
        Get the ids of all reference sequences in all alignment files.

        @param alignmentFiles: A C{list} of C{str} alignment file names.
        @return: A C{set} of C{str} reference ids as found in all passed
            alignment files.
        """
        # Get the names of all references in all alignment files.
        alignedReferences = set()
        for filename in alignmentFiles:
            with samfile(filename) as sam:
                for i in range(sam.nreferences):
                    alignedReferences.add(sam.get_reference_name(i))

        return alignedReferences

    def _readReferenceGenomes(self, referenceGenomeFiles):
        """
        Read reference genomes from files and check that any duplicates have
        identical sequences.

        @param referenceGenomeFiles: A C{list} of C{str} names of FASTA files
            containing reference genomes.
        @raise ValueError: If a reference genome is found in more than one file
            and the sequences are not identical.
        @return: A C{dict} keyed by C{str} sequence id with C{dark.Read}
            values holding reference genomes.
        """
        result = {}
        seen = {}
        for filename in referenceGenomeFiles:
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

        self.report(
            'Read %d reference genome%s:\n%s' % (
                len(result), s(len(result)),
                '\n'.join('  %s' % id_ for id_ in sorted(result))),
            requiredVerbosityLevel=2)

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

    def _removePreExistingTopLevelOutputDirFiles(self):
        """
        Remove all pre-existing files from the top-level output directory.
        """
        paths = list(map(str, chain(
            Path(self.outputDir).glob('result-summary.txt'))))

        if paths:
            self.report(
                '    Removing %d pre-existing output file%s from '
                'top-level output directory %s.' %
                (len(paths), s(len(paths)), self.outputDir),
                requiredVerbosityLevel=2)
            list(map(unlink, paths))

    def _removePreExistingAlignmentDirFiles(self, directory):
        """
        Remove all pre-existing files from the output directory for an
        alignment.

        @param directory: The C{str} directory to examine.
        """
        # This prevents us from doing a run that results in (say) 6
        # component files and then later doing a run that results in
        # only 5 components and erroneously thinking that
        # component-6-2.fasta etc. are from the most recent run.
        paths = list(map(str, chain(
            Path(directory).glob('*.stats'),
            Path(directory).glob('*.fasta'),
            Path(directory).glob('*.html'))))

        if paths:
            self.report(
                '    Removing %d pre-existing output file%s from %s '
                'directory.' % (len(paths), s(len(paths)), directory),
                requiredVerbosityLevel=2)
            list(map(unlink, paths))

    def _removePreExistingReferenceDirFiles(self, directory):
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
            self.report(
                '    Removing %d pre-existing output file%s from %s '
                'directory.' % (len(paths), s(len(paths)), directory),
                requiredVerbosityLevel=2)
            list(map(unlink, paths))

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

    def _plotCoverageAndSignificantLocations(
            self, referenceId, alignmentFile, readCountAtOffset, genomeLength,
            significantOffsets, outputDir):
        """
        Plot coverage and signifcant offsets.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFile: The C{str} name of an alignment file.
        @param readCountAtOffset: A C{list} of C{int} counts of the total
            number of reads at each genome offset (i.e., just the sum of the
            values in C{baseCountAtOffset})
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'coverage-and-significant-offsets.html')
        self.report('    Saving coverage and significant offset plot to',
                    filename)
        title = ('Coverage and significant offsets for aligment of %s in '
                 '%s' % (referenceId, alignmentFile))
        plotCoverageAndSignificantLocations(
            readCountAtOffset, genomeLength, significantOffsets, filename,
            title=title)

    def writeBcftoolsConsensus(self, referenceId, alignmentFile, outputDir):
        """
        Write a reference consensus using bcftools.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'reference-consensus-samtools.fasta')
        self.report('    Saving samtools reference consensus to', filename)

        referenceFilename = join(outputDir, 'reference.fasta')

        e = Executor()

        e.execute(
            'samtools mpileup -u -f %s %s 2>/dev/null | '
            'bcftools call -c | vcfutils.pl vcf2fq | '
            'filter-fasta.py --fastq --saveAs fasta --quiet '
            '--idLambda \'lambda _: "consensus-%s-samtools"\' > %s' %
            (referenceFilename, alignmentFile, referenceId, filename))

        if self.verbose > 1:
            for line in e.log:
                print('    ', line)

    def saveSignificantOffsets(self, significantOffsets, outputDir):
        """
        Save the significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'significant-offsets.txt')
        self.report('    Saving significant offsets to', filename)
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
        self.report('    Saving base nucleotide frequencies to', filename)

        genomeLengthWidth = int(log10(genomeLength)) + 1

        with open(filename, 'w') as fp:
            for offset in range(genomeLength):
                print('Location %*d: base counts %s' %
                      (genomeLengthWidth, offset + 1,
                       baseCountsToStr(baseCountAtOffset[offset])), file=fp)

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

        filename = join(outputDir, 'reference.fasta')
        self.report('    Saving reference to', filename)
        Reads([reference]).save(filename)

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

    def saveReferenceBaseFrequencyPlot(
            self, referenceId, genomeLength, significantOffsets,
            baseCountAtOffset, readCountAtOffset, outputDir):
        """
        Make a plot of the sorted base frequencies for the reference.

        @param referenceId: The C{str} id of the reference sequence.
        @param genomeLength: The C{int} length of the genome the reads were
            aligned to.
        @param significantOffsets: A C{set} of signifcant offsets.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param readCountAtOffset: A C{list} of C{int} counts of the total
            number of reads at each genome offset (i.e., just the sum of the
            values in C{baseCountAtOffset})
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, 'reference-base-frequencies.html')
        self.report('    Writing reference base frequency plot to', filename)
        plotBaseFrequencies(
            significantOffsets, baseCountAtOffset, readCountAtOffset, filename,
            title='%s (length %d)' % (referenceId, genomeLength),
            minReads=self.minReads, homogeneousCutoff=self.homogeneousCutoff,
            histogram=False, show=False)

    def saveReducedFasta(self, significantOffsets, outputDir):
        """
        Write out FASTA that contains reads with bases just at the
        significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        self.report('    Saving reduced FASTA')
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

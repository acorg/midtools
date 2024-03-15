from __future__ import print_function, division

import sys
from tempfile import mkdtemp
from os import unlink
from os.path import exists, join, basename
from os import mkdir
from math import log10
from pathlib import Path  # This is Python 3 only.
from itertools import chain
from collections import defaultdict

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Reads
from dark.sam import SAMFilter, PaddedSAM, samfile

from midtools.data import gatherData, findSignificantOffsets
from midtools.plotting import (
    plotBaseFrequencies,
    plotCoverageAndSignificantLocations,
    plotSAM,
)
from midtools.read import AlignedRead
from midtools.utils import baseCountsToStr, fastaIdentityTable, s, commas
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
    @param plotSAM: If C{True} save plots of where reads lie on each reference
        genome (can be slow).
    @param saveReducedFASTA: If C{True}, write out a FASTA file of the original
        input but with just the signifcant locations.
    @param verbose: The C{int}, verbosity level. Use C{0} for no output.
    """

    DEFAULT_HOMOGENEOUS_CUTOFF = 0.9
    DEFAULT_MIN_READS = 5

    def __init__(
        self,
        alignmentFiles,
        referenceGenomeFiles,
        referenceIds=None,
        outputDir=None,
        minReads=DEFAULT_MIN_READS,
        homogeneousCutoff=DEFAULT_HOMOGENEOUS_CUTOFF,
        plotSAM=False,
        saveReducedFASTA=False,
        verbose=0,
    ):
        self.alignmentFiles = alignmentFiles
        self.outputDir = outputDir
        self.minReads = minReads
        self.homogeneousCutoff = homogeneousCutoff
        self.plotSAM = plotSAM
        self.saveReducedFASTA = saveReducedFASTA
        self.verbose = verbose
        self.referenceGenomes = self._readReferenceGenomes(referenceGenomeFiles)

        # Make short reference ids from the reference genomes.
        self.shortReferenceId = dict(
            (id_, id_.split()[0]) for id_ in self.referenceGenomes
        )

        # Make short output file names from the given reference file names.
        self.shortAlignmentFilename = dict(
            (filename, basename(filename).rsplit(".", maxsplit=1)[0])
            for filename in alignmentFiles
        )

        alignedReferences = self._getAlignedReferences(alignmentFiles)
        self.referenceIds = self._getReferenceIds(alignedReferences, referenceIds)

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
                    "Alignments against the following reference id%s are not "
                    "present in any alignment file:\n%s"
                    % (
                        s(len(missing)),
                        "\n".join("  %s" % id_ for id_ in sorted(missing)),
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)

            missing = set(referenceIds) - set(self.referenceGenomes)
            if missing:
                print(
                    "Reference id%s %s not present in any reference genome "
                    "file." % (s(len(missing)), commas(missing)),
                    file=sys.stderr,
                )
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
                    "No analysis will be performed on reference%s %s "
                    "(found in SAM/BAM alignment file(s) headers) because no "
                    "corresponding reference genome was found."
                    % (s(len(missing)), commas(missing))
                )

            referenceIds = alignedReferences & set(self.referenceGenomes)

            if referenceIds:
                self.report(
                    "Examining %d reference%s: %s"
                    % (len(referenceIds), s(len(referenceIds)), commas(referenceIds))
                )
            else:
                print(
                    "Nothing to do! No genome could be found for any aligned "
                    "reference. Found reference%s: %s"
                    % (s(len(alignedReferences)), commas(alignedReferences)),
                    file=sys.stderr,
                )
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
            self.report("Analyzing alignment file", alignmentFile)
            alignmentOutputDir = self._setupAlignmentOutputDir(alignmentFile, outputDir)

            self._writeAlignmentFileSummary(alignmentFile, alignmentOutputDir)

            for referenceId in sorted(self.referenceIds):
                self.report("  Looking for reference", referenceId)

                referenceOutputDir = self._setupReferenceOutputDir(
                    referenceId, alignmentOutputDir
                )

                result = self.analyzeReferenceId(
                    referenceId, alignmentFile, referenceOutputDir
                )

                if result:
                    results[alignmentFile][referenceId] = result

            self._writeAlignmentHTMLSummary(results[alignmentFile], alignmentOutputDir)

        self._writeOverallResultSummary(results, self.outputDir)
        self._writeOverallResultSummarySummary(results, self.outputDir)

        return results

    def analyzeReferenceId(self, referenceId, alignmentFile, outputDir):
        """
        Analyze the given reference id in the given alignment file (if an
        alignment to the reference id is present).

        @param referenceId: The C{str} id of the reference sequence to analyze.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the output directory.
        @return: C{None} if C{referenceId} is not present in C{alignmentFile}
            or if no significant offsets are found. Else, a C{dict} with C{str}
            keys 'significantOffsets' (containing the signifcant offsets) and
            'consensusRead', the consensus sequence that best matches
            C{referenceId}.
        """
        raise NotImplementedError("Subclasses must implement this method")

    def _writeAlignmentFileSummary(self, alignmentFile, outputDir):
        """
        Write a summary of alignments.

        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{str} name of the output directory.
        """
        shortAlignmentFilename = self.shortAlignmentFilename[alignmentFile]
        filename = join(outputDir, shortAlignmentFilename + ".stats")
        self.report("  Writing alignment statistics to", filename)
        e = Executor()
        e.execute('sam-reference-read-counts.py "%s" > %s' % (alignmentFile, filename))
        if self.verbose > 1:
            for line in e.log:
                print("    ", line)

    def _writeAlignmentHTMLSummary(self, result, outputDir):
        """
        Write an HTML summary of the overall results.

        @param result: A C{dict} keyed by C{str} short reference name, and
           with values being C{dict}s with signifcant offsets and best
           consensus sequence for the corresponding reference in the alignment
           file.
        """
        referencesFilename = join(outputDir, "references.fasta")
        self.report("  Writing FASTA for mapped-to references to", referencesFilename)
        with open(referencesFilename, "w") as fp:
            for referenceId in sorted(result):
                print(
                    self.referenceGenomes[referenceId].toString("fasta"),
                    file=fp,
                    end="",
                )

        consensusesFilename = join(outputDir, "consensuses.fasta")
        self.report(
            "  Writing FASTA consensus for mapped-to references to", consensusesFilename
        )
        with open(consensusesFilename, "w") as fp:
            for referenceId in sorted(result):
                print(
                    result[referenceId]["consensusRead"].toString("fasta"),
                    file=fp,
                    end="",
                )

        htmlFilename = join(outputDir, "consensus-vs-reference.html")
        self.report("  Writing consensus vs reference identity table to", htmlFilename)
        fastaIdentityTable(
            consensusesFilename,
            htmlFilename,
            self.verbose,
            filename2=referencesFilename,
        )

        htmlFilename = join(outputDir, "consensus-vs-consensus.html")
        self.report("  Writing consensus vs consensus identity table to", htmlFilename)
        fastaIdentityTable(consensusesFilename, htmlFilename, self.verbose)

    def _writeOverallResultSummary(self, results, outputDir):
        """
        Write a summary of the overall results.

        @param results: A C{dict} of C{dicts}. Keyed by C{str} short alignment
           file name, then C{str} short reference name, and with values being
           C{dict}s with signifcant offsets and best consensus sequence for
           the corresponding reference in the alignment file.
        """
        filename = join(outputDir, "result-summary.txt")
        self.report("Writing overall result summary to", filename)
        with open(filename, "w") as fp:
            for alignmentFilename in sorted(results):
                print("Alignment file", alignmentFilename, file=fp)
                for referenceId in sorted(results[alignmentFilename]):
                    result = results[alignmentFilename][referenceId]
                    referenceRead = self.referenceGenomes[referenceId]
                    consensusRead = result["consensusRead"]
                    genomeLength = len(referenceRead)
                    significantOffsets = result["significantOffsets"]
                    print(
                        "\n  Reference %s (length %d)" % (referenceId, genomeLength),
                        file=fp,
                    )
                    print(
                        "    %d significant offsets found." % len(significantOffsets),
                        file=fp,
                    )

                    # Overall match.
                    match = compareDNAReads(referenceRead, consensusRead)
                    print("\n    Overall match of reference with consensus:", file=fp)
                    print(
                        matchToString(
                            match, referenceRead, consensusRead, indent="    "
                        ),
                        file=fp,
                    )

                    # Significant sites match.
                    match = compareDNAReads(
                        referenceRead, consensusRead, offsets=significantOffsets
                    )
                    print(
                        "\n    Match of reference with consensus at "
                        "%d SIGNIFICANT sites:" % len(significantOffsets),
                        file=fp,
                    )
                    print(
                        matchToString(
                            match,
                            referenceRead,
                            consensusRead,
                            indent="    ",
                            offsets=significantOffsets,
                        ),
                        file=fp,
                    )

                    # Non-significant sites match.
                    nonSignificantOffsets = set(range(genomeLength)) - set(
                        significantOffsets
                    )
                    match = compareDNAReads(
                        referenceRead, consensusRead, offsets=nonSignificantOffsets
                    )
                    print(
                        "\n    Match of reference with consensus at "
                        "%d NON-SIGNIFICANT sites:" % len(nonSignificantOffsets),
                        file=fp,
                    )
                    print(
                        matchToString(
                            match,
                            referenceRead,
                            consensusRead,
                            indent="    ",
                            offsets=nonSignificantOffsets,
                        ),
                        file=fp,
                    )

    def _writeOverallResultSummarySummary(self, results, outputDir):
        """
        Write a summary of the summary of the overall results.

        @param results: A C{dict} of C{dicts}. Keyed by C{str} short alignment
           file name, then C{str} short reference name, and with values being
           C{dict}s with signifcant offsets and best consensus sequence for
           the corresponding reference in the alignment file.
        """
        filename = join(outputDir, "result-summary-summary.txt")
        self.report("Writing overall result summary summary to", filename)

        bestFraction = 0.0
        bestAlignmentReference = []

        with open(filename, "w") as fp:
            for alignmentFilename in sorted(results):
                print(alignmentFilename, file=fp)
                resultSummary = []
                for referenceId in sorted(results[alignmentFilename]):
                    result = results[alignmentFilename][referenceId]
                    referenceRead = self.referenceGenomes[referenceId]
                    consensusRead = result["consensusRead"]
                    match = compareDNAReads(referenceRead, consensusRead)["match"]
                    matchCount = (
                        match["identicalMatchCount"] + match["ambiguousMatchCount"]
                    )
                    fraction = matchCount / len(referenceRead)

                    if fraction > bestFraction:
                        bestFraction = fraction
                        bestAlignmentReference = [(alignmentFilename, referenceId)]
                    elif fraction == bestFraction:
                        bestAlignmentReference.append((alignmentFilename, referenceId))

                    resultSummary.append(
                        (
                            fraction,
                            "  %s: %d/%d (%.2f%%)"
                            % (
                                referenceId,
                                matchCount,
                                len(referenceRead),
                                fraction * 100.0,
                            ),
                        )
                    )

                # Sort the result summary by decreasing nucleotide identity
                # fraction.
                resultSummary.sort(reverse=True)
                for fraction, summary in resultSummary:
                    print(summary, file=fp)

                print(file=fp)

            print(
                "Best match%s (%.2f%%):"
                % (
                    "" if len(bestAlignmentReference) == 1 else "es",
                    bestFraction * 100.0,
                ),
                file=fp,
            )
            for alignmentFilename, referenceId in bestAlignmentReference:
                print("  %s: %s" % (alignmentFilename, referenceId), file=fp)

    def initialReferenceIdAnalysis(self, referenceId, alignmentFile, outputDir):
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
                self.report("    Reference %s not in alignment file." % referenceId)
                return
            else:
                genomeLength = sam.lengths[tid]
                # Sanity check.
                assert genomeLength == len(self.referenceGenomes[referenceId])

        if self.plotSAM:
            filename = join(outputDir, "reads.html")
            self.report("    Saving reads alignment plot to %s" % filename)
            plotSAM(
                SAMFilter(alignmentFile, referenceIds={referenceId}),
                filename,
                title=referenceId,
                jitter=0.45,
            )

        alignedReads = []
        samFilter = SAMFilter(
            alignmentFile,
            referenceIds={referenceId},
            dropDuplicates=True,
            dropSupplementary=True,
            # dropSecondary=True,
            storeQueryIds=True,
        )
        paddedSAM = PaddedSAM(samFilter)
        for query in paddedSAM.queries(addAlignment=True):
            assert len(query) == genomeLength
            alignedReads.append(AlignedRead(query.id, query.sequence, query.alignment))

        # Sanity check that all aligned reads have different ids. This
        # should be the case because the padded SAM queries method adds /2,
        # /3 etc to queries that have more than one alignment.
        assert len(alignedReads) == len(set(read.id for read in alignedReads))

        readCountAtOffset, baseCountAtOffset, readsAtOffset = gatherData(
            genomeLength, alignedReads
        )

        significantOffsets = list(
            findSignificantOffsets(
                baseCountAtOffset,
                readCountAtOffset,
                self.minReads,
                self.homogeneousCutoff,
            )
        )

        self.report(
            "    %d alignment%s (of %d unique %s) read from %s"
            % (
                samFilter.alignmentCount,
                s(samFilter.alignmentCount),
                len(samFilter.queryIds),
                "query" if len(samFilter.queryIds) == 1 else "queries",
                alignmentFile,
            )
        )
        self.report(
            "    %d of which %s aligned to %s"
            % (
                len(alignedReads),
                "was" if len(alignedReads) == 1 else "were",
                referenceId,
            )
        )
        self.report("    Reference genome length %d" % genomeLength)
        self.report(
            "    Found %d significant location%s"
            % (len(significantOffsets), s(len(significantOffsets)))
        )

        self.saveBaseFrequencies(outputDir, genomeLength, baseCountAtOffset)

        if not significantOffsets:
            self.report("    No significant locations found.")
            return

        if self.saveReducedFASTA:
            self.saveReducedFasta(significantOffsets, outputDir)

        self._plotCoverageAndSignificantLocations(
            referenceId,
            alignmentFile,
            readCountAtOffset,
            genomeLength,
            significantOffsets,
            outputDir,
        )

        self.saveSignificantOffsets(significantOffsets, outputDir)

        for read in alignedReads:
            read.setSignificantOffsets(significantOffsets)

        self.saveReferenceBaseFrequencyPlot(
            referenceId,
            genomeLength,
            significantOffsets,
            baseCountAtOffset,
            readCountAtOffset,
            outputDir,
        )

        # Save the reference.
        filename = join(outputDir, "reference.fasta")
        self.report("    Saving reference to", filename)
        reference = self.referenceGenomes[referenceId]
        Reads([reference]).save(filename)

        # Extract a consensus according to bcftools.
        self.writeBcftoolsConsensus(referenceId, alignmentFile, outputDir)

        return (
            genomeLength,
            alignedReads,
            readCountAtOffset,
            baseCountAtOffset,
            readsAtOffset,
            significantOffsets,
            samFilter,
            paddedSAM,
        )

    def _setupOutputDir(self):
        """
        Set up the output directory and return its path.

        @return: The C{str} path of the output directory.
        """
        if self.outputDir:
            if exists(self.outputDir):
                self._removePreExistingTopLevelOutputDirFiles()
            else:
                mkdir(self.outputDir)
        else:
            self.outputDir = mkdtemp()
            print("Writing output files to %s" % self.outputDir)
        return self.outputDir

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
                            "Reference genome id %r was found in two files "
                            "(%r and %r) but with different sequences."
                            % (id_, seen[id_], filename)
                        )
                else:
                    seen[id_] = filename
                    result[id_] = read

        self.report(
            "Read %d reference genome%s:\n%s"
            % (
                len(result),
                s(len(result)),
                "\n".join("  %s" % id_ for id_ in sorted(result)),
            ),
            requiredVerbosityLevel=2,
        )

        return result

    def _removePreExistingTopLevelOutputDirFiles(self):
        """
        Remove all pre-existing files from the top-level output directory.
        """
        paths = list(map(str, chain(Path(self.outputDir).glob("result-summary.txt"))))

        if paths:
            self.report(
                "    Removing %d pre-existing output file%s from "
                "top-level output directory %s."
                % (len(paths), s(len(paths)), self.outputDir),
                requiredVerbosityLevel=2,
            )
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
        paths = list(
            map(
                str,
                chain(
                    Path(directory).glob("*.stats"),
                    Path(directory).glob("*.fasta"),
                    Path(directory).glob("*.html"),
                ),
            )
        )

        if paths:
            self.report(
                "    Removing %d pre-existing output file%s from %s "
                "directory." % (len(paths), s(len(paths)), directory),
                requiredVerbosityLevel=2,
            )
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
        paths = list(
            map(
                str,
                chain(
                    Path(directory).glob("*.fasta"),
                    Path(directory).glob("*.html"),
                    Path(directory).glob("*.txt"),
                ),
            )
        )

        if paths:
            self.report(
                "    Removing %d pre-existing output file%s from %s "
                "directory." % (len(paths), s(len(paths)), directory),
                requiredVerbosityLevel=2,
            )
            list(map(unlink, paths))

    def _plotCoverageAndSignificantLocations(
        self,
        referenceId,
        alignmentFile,
        readCountAtOffset,
        genomeLength,
        significantOffsets,
        outputDir,
    ):
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
        filename = join(outputDir, "coverage-and-significant-offsets.html")
        self.report("    Saving coverage and significant offset plot to", filename)
        title = "Coverage and significant offsets for alignment of %s in " "%s" % (
            referenceId,
            alignmentFile,
        )
        plotCoverageAndSignificantLocations(
            readCountAtOffset, genomeLength, significantOffsets, filename, title=title
        )

    def writeBcftoolsConsensus(self, referenceId, alignmentFile, outputDir):
        """
        Write a reference consensus using bcftools.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, "reference-consensus-samtools.fasta")
        self.report("    Saving samtools reference consensus to", filename)
        referenceFilename = join(outputDir, "reference.fasta")

        e = Executor()

        e.execute(
            "samtools mpileup -u -f %s %s 2>/dev/null | "
            "bcftools call -c | vcfutils.pl vcf2fq | "
            "filter-fasta.py --fastq --saveAs fasta --quiet "
            "--idLambda 'lambda _: \"consensus-%s-samtools\"' > %s"
            % (referenceFilename, alignmentFile, referenceId, filename)
        )

        if self.verbose > 1:
            for line in e.log:
                print("    ", line)

    def saveSignificantOffsets(self, significantOffsets, outputDir):
        """
        Save the significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        filename = join(outputDir, "significant-offsets.txt")
        self.report("    Saving significant offsets to", filename)
        with open(filename, "w") as fp:
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
        filename = join(outputDir, "base-frequencies.txt")
        self.report("    Saving base nucleotide frequencies to", filename)

        genomeLengthWidth = int(log10(genomeLength)) + 1

        with open(filename, "w") as fp:
            for offset in range(genomeLength):
                print(
                    "Location %*d: base counts %s"
                    % (
                        genomeLengthWidth,
                        offset + 1,
                        baseCountsToStr(baseCountAtOffset[offset]),
                    ),
                    file=fp,
                )

    def saveReferenceBaseFrequencyPlot(
        self,
        referenceId,
        genomeLength,
        significantOffsets,
        baseCountAtOffset,
        readCountAtOffset,
        outputDir,
    ):
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
        filename = join(outputDir, "reference-base-frequencies.html")
        self.report("    Writing reference base frequency plot to", filename)
        plotBaseFrequencies(
            significantOffsets,
            baseCountAtOffset,
            readCountAtOffset,
            filename,
            title="%s (length %d)" % (referenceId, genomeLength),
            minReads=self.minReads,
            homogeneousCutoff=self.homogeneousCutoff,
            histogram=False,
            show=False,
        )

    def saveReducedFasta(self, significantOffsets, outputDir):
        """
        Write out FASTA that contains reads with bases just at the
        significant offsets.

        @param significantOffsets: A C{set} of signifcant offsets.
        @param outputDir: A C{str} directory path.
        """
        self.report("    Saving reduced FASTA")
        print("    Saving reduced FASTA not implemented yet")
        return

        allGaps = "-" * len(significantOffsets)

        def unwanted(read):
            return None if read.sequence == allGaps else read

        FastaReads(self.fastaFile).filter(keepSites=significantOffsets).filter(
            modifier=unwanted
        ).save(join(outputDir, "reduced.fasta"))

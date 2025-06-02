from tempfile import mkdtemp
from os import unlink
from pathlib import Path
from itertools import chain
from collections import defaultdict
from typing import Optional, Callable

from dark.dna import compareDNAReads
from dark.process import Executor
from dark.utils import pct

try:
    # This relies on the private Cambridge/Charite HBV repo.
    from pyhbv.genotype import getGenotype
except ImportError:
    haveHBV = False
    getGenotype = None
else:
    haveHBV = True

from midtools.match import matchToString
from midtools.reference import (
    getReferenceIds,
    getReferenceLength,
    Reference,
    readReferenceGenomes,
)
from midtools.utils import fastaIdentityTable, s, quoted


class ReadAnalysis:
    """
    Perform a read alignment analysis for multiple infection detection.

    @param sampleName: The C{str} name of the sample whose reads are being
        analysed.
    @param alignmentFiles: A C{list} of C{str} names of SAM/BAM alignment
        files. These files should have mappings of reads to all references
        (as produced using the --all option when running bowtie2).
    @param referenceGenomeFiles: A C{list} of C{str} names of FASTA files
        containing reference genomes.
    @param outputDir: The C{str} directory to save result files to.
    @param referenceIds: The C{str} sequence ids whose alignment should be
        analyzed. All ids must be present in the C{referenceGenomes} files.
        One of the SAM/BAM files given using C{alignmentFiles} should have an
        alignment against the given argument. If omitted, all references that
        are aligned to in the given BAM/SAM files will be analyzed.
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
        sampleName: str,
        alignmentFiles: list[str],
        referenceGenomeFiles: list[str],
        outputDir: str,
        referenceIds: Optional[list[str]] = None,
        minReads: int = DEFAULT_MIN_READS,
        homogeneousCutoff: float = DEFAULT_HOMOGENEOUS_CUTOFF,
        plotSAM: bool = False,
        saveReducedFASTA: bool = False,
        verbose: int = 0,
    ) -> None:
        self.sampleName = sampleName
        self.alignmentFiles = [Path(f) for f in alignmentFiles]
        self.outputDir = Path(outputDir)
        self.minReads = minReads
        self.homogeneousCutoff = homogeneousCutoff
        self.plotSAM = plotSAM
        self.saveReducedFASTA = saveReducedFASTA
        self.verbose = verbose
        self.referenceGenomes = readReferenceGenomes(referenceGenomeFiles)

        # Make short reference ids from the reference genomes.
        self.shortReferenceId = dict(
            (id_, id_.split()[0]) for id_ in self.referenceGenomes
        )

        # Make short output file names from the given reference file names.
        self.shortAlignmentFilename = dict(
            (filename, filename.name.rsplit(".", maxsplit=1)[0])
            for filename in self.alignmentFiles
        )

        self.referenceIds = getReferenceIds(self, referenceIds)

    def report(self, *args, requiredVerbosityLevel: int = 1) -> None:
        """
        Print a status message, if our verbose setting is high enough.

        @param args: The arguments to print.
        @param requiredVerbosityLevel: The minimum C{int} verbosity
            level required.
        """
        if self.verbose >= requiredVerbosityLevel:
            print(*args)

    def run(
        self, analysisFunc: Callable[[Reference], None]
    ) -> dict[Path, dict[str, Reference]]:
        """
        Perform a read analysis for all reference sequences.

        @param analysisFunc: A function that takes a C{Reference} instance and
            performs specialized analysis on it.
        """
        outputDir = self._setupOutputDir()
        results: dict[str, dict[str, Reference]] = defaultdict(
            lambda: defaultdict(dict)
        )

        for alignmentFile in self.alignmentFiles:
            self.report(f"Analyzing alignment file {quoted(alignmentFile)}")
            alignmentOutputDir = self._setupAlignmentOutputDir(alignmentFile, outputDir)

            self._writeAlignmentFileSummary(alignmentFile, alignmentOutputDir)

            for count, referenceId in enumerate(sorted(self.referenceIds), start=1):
                self.report(
                    f"  Looking for reference {referenceId!r} "
                    f"({count}/{len(self.referenceIds)})."
                )

                referenceLength = getReferenceLength(referenceId, alignmentFile)

                if referenceLength is None:
                    self.report(f"    Reference {referenceId!r} not in alignment file.")
                    continue

                referenceOutputDir = self._setupReferenceOutputDir(
                    referenceId, alignmentOutputDir
                )

                reference = Reference(
                    self,
                    referenceId,
                    self.referenceGenomes[referenceId],
                    alignmentFile,
                    referenceOutputDir,
                )

                analysisFunc(reference)

                results[alignmentFile][referenceId] = reference

            self._writeAlignmentHTMLSummary(results[alignmentFile], alignmentOutputDir)

        self._writeOverallResultSummary(results, outputDir)
        self._writeOverallResultSummarySummary(results, outputDir)

        return results

    def analyzeReference(self, reference: Reference):
        """
        Analyze the given reference id in the given alignment file (if an
        alignment to the reference id is present).

        @param referenceId: The C{str} id of the reference sequence to analyze.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{Path} to the output directory.
        @return: C{None} if C{referenceId} is not present in C{alignmentFile}
            or if no significant offsets are found. Else, a C{dict} with C{str}
            keys 'significantOffsets' (containing the signifcant offsets) and
            'consensusRead', the consensus sequence that best matches
            C{referenceId}.
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def _writeAlignmentFileSummary(self, alignmentFile: Path, outputDir: Path) -> None:
        """
        Write a summary of alignments.

        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{Path} to the output directory.
        """
        shortAlignmentFilename = self.shortAlignmentFilename[alignmentFile]
        filename = outputDir / (shortAlignmentFilename + ".stats")
        self.report("  Writing alignment statistics to", filename)
        e = Executor()
        e.execute(
            f"sam-reference-read-counts.py {quoted(alignmentFile)} > {quoted(filename)}"
        )
        if self.verbose > 1:
            for line in e.log:
                print("    ", line)

    def _writeAlignmentHTMLSummary(
        self, result: dict[str, dict[str, Reference]], outputDir: Path
    ) -> None:
        """
        Write an HTML summary of the overall results.

        @param result: A C{dict} keyed by C{str} short reference name, and
           with values being C{dict}s with signifcant offsets and best
           consensus sequence for the corresponding reference in the alignment
           file.
        @param outputDir: The C{Path} to the output directory.
        """
        referencesFilename = outputDir / "references.fasta"
        self.report("  Writing FASTA for mapped-to references to", referencesFilename)
        with open(referencesFilename, "w") as fp:
            for referenceId in sorted(result):
                print(
                    self.referenceGenomes[referenceId].toString("fasta"),
                    file=fp,
                    end="",
                )

        consensusesFilename = outputDir / "consensuses.fasta"
        self.report(
            "  Writing FASTA consensus for mapped-to references to", consensusesFilename
        )
        with open(consensusesFilename, "w") as fp:
            for referenceId in sorted(result):
                print(
                    result[referenceId].consensusRead.toString("fasta"),
                    file=fp,
                    end="",
                )

        htmlFilename = outputDir / "consensus-vs-reference.html"
        self.report("  Writing consensus vs reference identity table to", htmlFilename)
        fastaIdentityTable(
            consensusesFilename,
            htmlFilename,
            self.verbose,
            filename2=referencesFilename,
        )

        htmlFilename = outputDir / "consensus-vs-consensus.html"
        self.report("  Writing consensus vs consensus identity table to", htmlFilename)
        fastaIdentityTable(consensusesFilename, htmlFilename, self.verbose)

    def _writeOverallResultSummary(
        self, results: dict[Path, dict[str, Reference]], outputDir: Path
    ) -> None:
        """
        Write a summary of the overall results.

        @param results: A C{dict} of C{dicts}. Keyed by C{str} short alignment
           file name, then C{str} short reference name, and with values being
           C{dict}s with signifcant offsets and best consensus sequence for
           the corresponding reference in the alignment file.
        @param outputDir: The C{Path} to the output directory.
        """
        filename = outputDir / "result-summary.txt"
        self.report("Writing overall result summary to", filename)
        with open(filename, "w") as fp:
            for alignmentFilename in sorted(results):
                print(f"Alignment file {quoted(alignmentFilename)}", file=fp)
                for referenceId in sorted(results[alignmentFilename]):
                    result = results[alignmentFilename][referenceId]
                    referenceRead = self.referenceGenomes[referenceId]
                    consensusRead = result.consensusRead
                    genomeLength = len(referenceRead)
                    significantOffsets = result.significantOffsets
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

    def _writeOverallResultSummarySummary(
        self, results: dict[Path, dict[str, Reference]], outputDir: Path
    ) -> None:
        """
        Write a summary of the summary of the overall results.

        @param results: A C{dict} of C{dicts}. Keyed by C{str} short alignment
           file name, then C{str} short reference name, and with values being
           C{dict}s with signifcant offsets and best consensus sequence for
           the corresponding reference in the alignment file.
        @param outputDir: The C{Path} to the output directory.
        """
        filename = outputDir / "result-summary-summary.txt"
        self.report("Writing overall result summary summary to", filename)

        with open(filename, "w") as fp:
            for alignmentFilename in sorted(results):
                print(f"{str(alignmentFilename)}", file=fp)
                resultSummary = []
                for referenceId in sorted(results[alignmentFilename]):
                    result = results[alignmentFilename][referenceId]
                    referenceRead = self.referenceGenomes[referenceId]
                    consensusRead = result.consensusRead
                    match = compareDNAReads(referenceRead, consensusRead)["match"]
                    strictCount = match["identicalMatchCount"]
                    nonStrictCount = strictCount + match["ambiguousMatchCount"]
                    if haveHBV:
                        assert getGenotype
                        summary = "  %-20s\t%-10s\tStrict: %s\tNon-strict: %s" % (
                            referenceId,
                            getGenotype(referenceId),
                            pct(strictCount, len(referenceRead)),
                            pct(nonStrictCount, len(referenceRead)),
                        )
                    else:
                        summary = "  %-20s\tStrict: %s\tNon-strict: %s" % (
                            referenceId,
                            pct(strictCount, len(referenceRead)),
                            pct(nonStrictCount, len(referenceRead)),
                        )
                    resultSummary.append((strictCount, nonStrictCount, summary))

                # Sort the result summary by decreasing nucleotide identity
                # fraction. The strict count is first in the tuples we're sorting.
                resultSummary.sort(reverse=True)
                for _, _, summary in resultSummary:
                    print(summary, file=fp)

                print(file=fp)

    def _setupOutputDir(self) -> Path:
        """
        Set up the output directory and return its path.

        @return: The C{str} path of the output directory.
        """
        if self.outputDir:
            if self.outputDir.exists():
                self._removePreExistingTopLevelOutputDirFiles()
            else:
                self.outputDir.mkdir()
        else:
            self.outputDir = Path(mkdtemp())
            print("Writing output files to %s" % self.outputDir)
        return self.outputDir

    def _alignmentOutputDir(self, alignmentFile: Path, outputDir: Path) -> Path:
        """
        Get output directory for a given alignment file.

        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{Path} to the top-level output directory.
        @return: The C{str} output directory name.
        """
        return outputDir / self.shortAlignmentFilename[alignmentFile]

    def _setupAlignmentOutputDir(self, alignmentFile: Path, outputDir: Path) -> Path:
        """
        Set up the output directory for a given alignment file.

        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{Path} to the top-level output directory.
        @return: The C{str} output directory name.
        """
        directory = self._alignmentOutputDir(alignmentFile, outputDir)
        if directory.exists():
            self._removePreExistingAlignmentDirFiles(directory)
        else:
            directory.mkdir()

        return directory

    def _setupReferenceOutputDir(self, referenceId, outputDir):
        """
        Set up the output directory for a given alignment file and reference.

        @param referenceId: The C{str} id of the reference sequence.
        @param outputDir: The C{Path} to the top-level output directory.
        @return: The C{str} output directory name.
        """
        # Make short versions of the reference id and filename for a
        # per-alignment-file per-reference-sequence output directory.

        shortReferenceId = self.shortReferenceId[referenceId]
        directory = outputDir / shortReferenceId
        if directory.exists():
            self._removePreExistingReferenceDirFiles(directory)
        else:
            directory.mkdir()

        return directory

    def _removePreExistingTopLevelOutputDirFiles(self) -> None:
        """
        Remove all pre-existing files from the top-level output directory.
        """
        paths = list(map(str, chain(self.outputDir.glob("result-summary.txt"))))

        if paths:
            self.report(
                "    Removing %d pre-existing output file%s from "
                "top-level output directory %s."
                % (len(paths), s(len(paths)), str(self.outputDir)),
                requiredVerbosityLevel=2,
            )
            list(map(unlink, paths))

    def _removePreExistingAlignmentDirFiles(self, directory: Path) -> None:
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

    def _removePreExistingReferenceDirFiles(self, directory: Path) -> None:
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

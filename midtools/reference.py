import sys
from math import log10

from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Reads
from dark.sam import SAMFilter, PaddedSAM, samfile
from dark.utils import pct

from midtools.offsets import analyzeOffets, findSignificantOffsets
from midtools.plotting import (
    plotBaseFrequencies,
    plotCoverageAndSignificantLocations,
    plotSAM,
)
from midtools.read import AlignedRead
from midtools.utils import baseCountsToStr, commas, quoted, s


def _getAlignedReferenceIds(alignmentFiles):
    """
    Get the ids of all reference sequences in all alignment files.

    @param alignmentFiles: A C{list} of C{Path} alignment file names.
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


def getReferenceIds(readAnalysis, referenceIds):
    """
    Figure out which reference ids we can process.

    @param readAnalysis: A C{ReadAnalysis} instance.
    @param referenceIds: A C{list} of C{str} reference ids for which
        processing has specifically been requested, or C{None}.
    @return: A C{set} of C{str} reference ids to process.
    """

    alignedReferences = _getAlignedReferenceIds(readAnalysis.alignmentFiles)
    if referenceIds:
        # Specific reference ids were given. Check that each appears in some
        # alignment file and that we have a genome for each. Any error here causes a
        # message to stderr and exit.
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

        missing = set(referenceIds) - set(readAnalysis.referenceGenomes)
        if missing:
            print(
                "Reference id%s %s not present in any reference genome "
                "file." % (s(len(missing)), commas(missing)),
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        # We weren't told which reference ids to specifically examine the alignments
        # of, so examine all available references mentioned in any alignment file
        # and that we also have a genome for. Mention any references from alignment
        # files that we can't process due to lack of a genome.
        missing = alignedReferences - set(readAnalysis.referenceGenomes)
        if missing:
            readAnalysis.report(
                "No analysis will be performed on reference%s %s "
                "(found in SAM/BAM alignment file(s) headers) because no "
                "corresponding reference genome was found."
                % (s(len(missing)), commas(missing))
            )

        referenceIds = alignedReferences & set(readAnalysis.referenceGenomes)

        if referenceIds:
            readAnalysis.report(
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


def getReferenceLength(referenceId, alignmentFile):
    """
    Look for a reference name in a BAM file and return its length.

    @param referenceId: The C{str} id of the reference sequence.
    @param alignmentFile: The C{Path} to the BAM file supposedly containing matches to
        this reference.
    @return: An C{int} reference genome length.
    """

    with samfile(alignmentFile) as sam:
        tid = sam.get_tid(referenceId)
        if tid == -1:
            return
        else:
            return sam.lengths[tid]


def readReferenceGenomes(referenceGenomeFiles):
    """
    Read reference genomes from files and check that any duplicates have
    identical sequences.

    @param referenceGenomeFiles: A C{list} of C{Path} instances of FASTA
        files containing reference genomes.
    @raise ValueError: If a reference genome is found in more than one file
        and the sequences are not identical.
    @return: A C{dict} keyed by shortened C{str} sequence id with C{dark.Read}
        values holding reference genomes.
    """
    result = {}
    seen = {}
    for filename in referenceGenomeFiles:
        for read in FastaReads(filename):
            shortId = read.id.split()[0]
            if shortId in seen:
                if result[shortId].sequence != read.sequence:
                    raise ValueError(
                        "Reference genome id %r was found in two files "
                        "(%r and %r) but with different sequences."
                        % (shortId, seen[shortId], filename)
                    )
            else:
                seen[shortId] = filename
                result[shortId] = read

    return result


class Reference:
    """
    Hold information about an analyzed reference.

    @param readAnalysis: The C{ReadAnalysis} instance that that this reference is
        being analyzed in the context of.
    @param referenceId: The C{str} id of the reference sequence to analyze. This will
        typically be the short form (up to the first space) of C{read.id} (which
        contains the full reference id, as read from a FASTA file).
    @param read: A C{dark.reads.DNARead} instance.
    @param referenceLength: The C{int} length of the reference sequence.
    @param alignmentFile: The C{Path} to the BAM file supposedly containing matches to
        this reference.
    @param outputDir: The C{Path} to the output directory.
    """

    def __init__(self, readAnalysis, referenceId, read, alignmentFile, outputDir):
        self.readAnalysis = readAnalysis
        self.id = referenceId
        self.read = read
        self.alignmentFile = alignmentFile
        self.outputDir = outputDir
        self._initialAnalysis()

    def _initialAnalysis(self):
        """
        Analyze the given reference id in the given alignment file (if an
        alignment to the reference id is present).

        @return: C{None} if C{referenceId} is not present in C{alignmentFile}
            or if no significant offsets are found. Else, a C{dict} containing
            the signifcant offsets and the consensus sequence that best matches
            C{referenceId}.
        """

        report = self.readAnalysis.report

        if self.readAnalysis.plotSAM:
            filename = self.outputDir / "reads.html"
            report(f"    Saving reads alignment plot to {str(filename)!r}")
            plotSAM(
                SAMFilter(
                    self.alignmentFile,
                    referenceIds={self.id},
                ),
                filename,
                title=self.id,
                jitter=0.45,
            )

        samFilter = SAMFilter(
            self.alignmentFile,
            referenceIds={self.id},
            dropDuplicates=True,
            dropSupplementary=True,
            storeQueryIds=True,
        )
        self.paddedSAM = PaddedSAM(samFilter)

        self.alignedReads = []
        for query in self.paddedSAM.queries(addAlignment=True):
            assert len(query) == len(self.read)
            self.alignedReads.append(
                AlignedRead(query.id, query.sequence, query.alignment)
            )

        # Sanity check that all aligned reads have different ids. This
        # should be the case because the padded SAM queries method adds /2,
        # /3 etc to queries that have more than one alignment.
        assert len(self.alignedReads) == len(set(read.id for read in self.alignedReads))

        (
            self.readCountAtOffset,
            self.baseCountAtOffset,
            self.readsAtOffset,
        ) = analyzeOffets(len(self.read), self.alignedReads)

        self.significantOffsets = list(
            findSignificantOffsets(
                self.baseCountAtOffset,
                self.readCountAtOffset,
                self.readAnalysis.minReads,
                self.readAnalysis.homogeneousCutoff,
            )
        )

        report(
            "    %d alignment%s (of %d unique %s) read from %r."
            % (
                samFilter.alignmentCount,
                s(samFilter.alignmentCount),
                len(samFilter.queryIds),
                "query" if len(samFilter.queryIds) == 1 else "queries",
                str(self.alignmentFile),
            )
        )
        report(
            "    %d of which %s aligned to %s"
            % (
                len(self.alignedReads),
                "was" if len(self.alignedReads) == 1 else "were",
                self.id,
            )
        )
        report("    Reference genome length %d" % len(self.read))

        coveredOffsetCount = sum(bool(counts) for counts in self.baseCountAtOffset)
        report(f"    Covered offsets {pct(coveredOffsetCount, len(self.read))}")

        report(
            "    Found %d significant location%s"
            % (len(self.significantOffsets), s(len(self.significantOffsets)))
        )

        self.saveBaseFrequencies()

        if not self.significantOffsets:
            report("    No significant locations found.")
            return

        if self.readAnalysis.saveReducedFASTA:
            self.saveReducedFasta()

        self._plotCoverageAndSignificantLocations()

        self.saveSignificantOffsets()

        for read in self.alignedReads:
            read.setSignificantOffsets(self.significantOffsets)

        self.saveReferenceBaseFrequencyPlot()

        # Save the reference.
        filename = self.outputDir / "reference.fasta"
        report("    Saving reference to", filename)
        Reads([self.read]).save(filename)

        # Extract a consensus according to bcftools.
        self.writeBcftoolsConsensus()

    def _plotCoverageAndSignificantLocations(self):
        """
        Plot coverage and signifcant offsets.
        """
        filename = self.outputDir / "coverage-and-significant-offsets.html"
        self.readAnalysis.report(
            "    Saving coverage and significant offset plot to", str(filename)
        )
        title = "Coverage and significant offsets for alignment of %s in %s" % (
            self.id,
            str(self.alignmentFile),
        )
        plotCoverageAndSignificantLocations(
            self.readCountAtOffset,
            len(self.read),
            self.significantOffsets,
            filename,
            title=title,
        )

    def writeBcftoolsConsensus(self):
        """
        Write a reference consensus using bcftools.

        @param referenceId: The C{str} id of the reference sequence.
        @param alignmentFile: The C{str} name of an alignment file.
        @param outputDir: The C{Path} to the output directory.
        """
        filename = self.outputDir / "reference-consensus-samtools.fasta"
        self.readAnalysis.report("    Saving samtools reference consensus to", filename)
        referenceFilename = self.outputDir / "reference.fasta"

        e = Executor()

        e.execute(
            f"make-consensus.py --reference {quoted(referenceFilename)} "
            f"--bam {quoted(self.alignmentFile)} | "
            f"filter-fasta.py --quiet "
            f"--idLambda 'lambda _: \"consensus-{self.id}-samtools\"' "
            f"> {quoted(filename)}"
        )

        if self.readAnalysis.verbose > 1:
            for line in e.log:
                print("    ", line)

    def saveSignificantOffsets(self):
        """
        Save the significant offsets.
        """
        filename = self.outputDir / "significant-offsets.txt"
        self.readAnalysis.report("    Saving significant offsets to", filename)
        with open(filename, "w") as fp:
            for offset in self.significantOffsets:
                print(offset, file=fp)

    def saveBaseFrequencies(self):
        """
        Save the base nucleotide frequencies.
        """
        filename = self.outputDir / "base-frequencies.txt"
        self.readAnalysis.report("    Saving base nucleotide frequencies to", filename)

        referenceLengthWidth = int(log10(len(self.read))) + 1

        with open(filename, "w") as fp:
            for offset in range(len(self.read)):
                print(
                    "Location %*d: base counts %s"
                    % (
                        referenceLengthWidth,
                        offset + 1,
                        baseCountsToStr(self.baseCountAtOffset[offset]),
                    ),
                    file=fp,
                )

    def saveReferenceBaseFrequencyPlot(self):
        """
        Make a plot of the sorted base frequencies for the reference.
        """
        filename = self.outputDir / "reference-base-frequencies.html"
        self.readAnalysis.report(
            "    Writing reference base frequency plot to", filename
        )
        plotBaseFrequencies(
            len(self.read),
            self.significantOffsets,
            self.baseCountAtOffset,
            self.readCountAtOffset,
            filename,
            title="%s (length %d)" % (self.id, len(self.read)),
            minReads=self.readAnalysis.minReads,
            homogeneousCutoff=self.readAnalysis.homogeneousCutoff,
            histogram=False,
            show=False,
            yRange=None,
        )

    def saveReducedFasta(self):
        """
        Write out FASTA that contains reads with bases just at the
        significant offsets.
        """
        self.readAnalysis.report("    Saving reduced FASTA")
        print("    Saving reduced FASTA not implemented yet")
        return

        allGaps = "-" * len(self.significantOffsets)

        def unwanted(read):
            return None if read.sequence == allGaps else read

        FastaReads(self.fastaFile).filter(keepSites=self.significantOffsets).filter(
            modifier=unwanted
        ).save(self.outputDir / "reduced.fasta")

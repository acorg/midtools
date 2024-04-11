from collections import defaultdict, Counter

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.reads import Read, Reads
from dark.utils import pct

from midtools.analysis import ReadAnalysis
from midtools.clusters import ReadCluster, ReadClusters
from midtools.match import matchToString
from midtools.offsets import analyzeOffets, findSignificantOffsets
from midtools.plotting import plotBaseFrequencies, plotConsistentComponents
from midtools.utils import (
    alignmentQuality,
    baseCountsToStr,
    commas,
    commonest,
    fastaIdentityTable,
    nucleotidesToStr,
    quoted,
    s,
)


def connectedComponentsByOffset(significantReads, maxClusterDist, fp, verbose=0):
    """
    Yield sets of reads that are connected according to what significant
    offsets they cover (the nucleotides at those offsets are irrelevant at
    this point).

    @param significantReads: A C{set} of C{AlignedRead} instances, all of
        which cover at least one significant offset.
    @param maxClusterDist: A C{float} indicating the distance beyond which
        clustering should be aborted.
    @param fp: A file-like object to write information to.
    @param verbose: The C{int} verbosity level. Use C{0} for no output.
    @return: A generator that yields C{Component} instances.
    """
    while significantReads:
        significantRead = sorted(significantReads)[0]
        significantReads.remove(significantRead)
        component = {significantRead}
        offsets = set(significantRead.significantOffsets)
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
        yield Component(component, offsets, maxClusterDist, fp, verbose)


class ConsistentComponent:
    """
    Hold information about a set of reads that share significant offsets
    and which (largely) agree on the nucleotides present at those offsets.

    @param readCluster: A C{ReadCluster} instance. If C{None}, a new C{ReadCluster}
        will be created.
    """

    def __init__(self, readCluster=None):
        self.readCluster = readCluster or ReadCluster()

    def __len__(self):
        return len(self.readCluster.reads)

    @property
    def reads(self):
        return self.readCluster.reads

    @property
    def nucleotides(self):
        return self.readCluster.nucleotides

    def add(self, read):
        """
        Add a read to this cluster component.

        @param read: An C{alignedRead} instance.
        """
        return self.readCluster.add(read)

    @property
    def offsets(self):
        """
        Get the set of significant offsets covered by the reads in this component.

        @return: A C{set} of C{int} site offsets.
        """
        # This is a property to match the property in the Component class.
        return self.readCluster.offsets

    def update(self, reads):
        """
        Add reads to this cluster component.

        @param reads: An iterable of C{alignedRead} instances.
        """
        for read in reads:
            self.add(read)

    def saveFasta(self, fp):
        """
        Save all reads as FASTA.

        @param fp: A file pointer to write to.
        """
        for read in sorted(self.readCluster.reads):
            print(read.toString("fasta"), end="", file=fp)

    def savePaddedFasta(self, fp):
        """
        Save all reads as FASTA, padded with gaps to preserve alignment.

        @param fp: A file pointer to write to.
        """
        for read in sorted(self.readCluster.reads):
            print(read.toPaddedString(), end="", file=fp)

    def consensusBase(self, offset, referenceSequence, infoFp):
        """
        Get a consensus sequence.

        @param componentOffsets: The C{set} of offsets in this component. This is *not*
            the same as the offsets in this consistent component because this consistent
            component may not have reads for all offsets.
        @param referenceSequence: The C{str} reference sequence.
        @param infoFp: A file pointer to write draw (and other) info to.
        @return: A C{str} consensus sequence.
        """
        if offset in self.readCluster.nucleotides:
            referenceBase = referenceSequence[offset]
            base = commonest(
                self.readCluster.nucleotides[offset],
                referenceBase,
                drawFp=infoFp,
                drawMessage=(
                    f"WARNING: consensus draw at offset {offset}" + " %(baseCounts)s."
                ),
            )
        else:
            base = "-"

        return base

    def consensusSequence(self, componentOffsets, referenceSequence, infoFp):
        """
        Get a consensus sequence.

        @param componentOffsets: The C{set} of offsets in this component. This is *not*
            the same as the offsets in this consistent component because this consistent
            component may not have reads for all offsets.
        @param referenceSequence: The C{str} reference sequence.
        @param infoFp: A file pointer to write draw (and other) info to.
        @return: A C{str} consensus sequence.
        """
        return "".join(
            self.consensusBase(offset, referenceSequence, infoFp)
            for offset in sorted(componentOffsets)
        )

    def saveConsensus(
        self, count, componentOffsets, referenceSequence, consensusFp, infoFp
    ):
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
        @param infoFp: A file pointer to write draw (and other) info to.
        """
        sequence = self.consensusSequence(componentOffsets, referenceSequence, infoFp)
        id_ = "consistent-component-%d-consensus (based on %d reads)" % (
            count,
            len(self.readCluster.reads),
        )
        print(Read(id_, sequence).toString("fasta"), file=consensusFp, end="")

    def summarize(self, fp, count, componentOffsets, referenceSequence):
        """
        Write out a summary of this consistent component.

        @param fp: The file pointer to write to.
        @param count: The C{int} number of this component.
        @param componentOffsets: The C{set} of offsets in this component.
        @param referenceSequence: The C{str} reference sequence.
        """
        plural = s(len(self.readCluster.reads))
        print(
            "    Component %d: %d read%s, covering %d offset%s"
            % (
                count,
                len(self.readCluster.reads),
                plural,
                len(self.readCluster.nucleotides),
                s(len(self.readCluster.nucleotides)),
            ),
            file=fp,
        )
        print("    Nucleotide counts for each offset:", file=fp)
        print(nucleotidesToStr(self.readCluster.nucleotides, "      "), file=fp)
        print(
            "    Consensus sequence: %s"
            % self.consensusSequence(componentOffsets, referenceSequence, fp),
            file=fp,
        )
        print("    Read%s:" % plural, file=fp)
        for read in sorted(self.readCluster.reads):
            print("     ", read, file=fp)


class Component:
    """
    Hold information about a set of reads that share significant offsets,
    regardless of the nucleotides present at those offsets.

    Provide a method to create a list of subsets of these reads
    (ConsistentComponent instances) that are (largely) consistent in the
    nucleotides at their shared offsets.

    @param reads: A C{list} of C{AlignedSegment} instances.
    @param offsets: A C{set} of significant offsets covered by C{reads}.
    @param maxClusterDist: The C{float} inter-cluster distance above which
        clusters are too different to merge.
    @param fp: A file-like object to write information to.
    @param verbose: The C{int} verbosity level. Use C{0} for no output.
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
        nSelfReads = len(self)
        nCcReads = sum(map(len, self.consistentComponents))
        assert nSelfReads == nCcReads, "%d != %d" % (nSelfReads, nCcReads)
        offsets = set()
        for read in self.reads:
            offsets.update(read.significantOffsets)
        assert offsets == self.offsets

    def summarize(self, fp, count, referenceSequence):
        """
        Write out a summary of this component.

        @param fp: The file pointer to write to.
        @param count: The C{int} number of this component.
        @param referenceSequence: The C{str} reference sequence.
        """
        ccLengths = ", ".join(
            str(length) for length in map(len, self.consistentComponents)
        )
        print(
            "component %d: %d reads, covering %d offsets, split into %d "
            "clusters of lengths %s."
            % (
                count,
                len(self),
                len(self.offsets),
                len(self.consistentComponents),
                ccLengths,
            ),
            file=fp,
        )
        print("  offsets:", commas(self.offsets), file=fp)
        for read in sorted(self.reads):
            print("  ", read, file=fp)

        for i, cc in enumerate(self.consistentComponents, start=1):
            print(file=fp)
            cc.summarize(fp, i, self.offsets, referenceSequence)

    def saveFasta(self, outputDir, count, verbose):
        """
        @param outputDir: The C{Path} to the output directory.
        @param count: The C{int} number of this component.
        @param verbose: The C{int} verbosity level. Use C{0} for no output.
        """
        for i, cc in enumerate(self.consistentComponents, start=1):
            filename = outputDir / f"component-{count}-{i}.fasta"
            if verbose > 1:
                print(f"      Saving component {count} {i} FASTA to {str(filename)!r}")
            with open(filename, "w") as fp:
                cc.saveFasta(fp)

            filename = outputDir / ("component-%d-%d-padded.fasta" % (count, i))
            if verbose > 1:
                print(
                    "      Saving component %d %d padded FASTA to" % (count, i),
                    filename,
                )
            with open(filename, "w") as fp:
                cc.savePaddedFasta(fp)

    def findClusters(self, fp, verbose):
        """
        Find clusters of reads up to the maximum cluster distance threshold in
        self.maxClusterDist according to what nucleotides they have at their
        significant offsets.

        @param fp: A file-like object to write information to.
        @param verbose: The C{int} verbosity level. Use C{0} for no output.
        @return: A generator that yields C{ConsistentComponent} instances.
        """
        readClusters = ReadClusters()

        # Using 'sorted' here makes the clustering result deterministic.
        for read in sorted(self.reads):
            readClusters.add(read)

        for readCluster in readClusters.analyze(self.maxClusterDist, fp):
            if verbose:
                print(
                    "      Found cluster with %d reads, covering %d offsets"
                    % (len(readCluster.reads), len(readCluster.nucleotides)),
                    file=fp,
                )
            yield ConsistentComponent(readCluster)

    def saveConsensuses(self, outputDir, count, referenceSequence, verbose):
        """
        Write a consensus sequence for each consistent component of this component.

        @param outputDir: The C{Path} to the output directory.
        @param count: The C{int} number of this component.
        @param referenceSequence: The C{str} reference sequence.
        @param verbose: The C{int} verbosity level. Use C{0} for no output.
        """
        consensusFilename = outputDir / ("component-%d-consensuses.fasta" % count)
        infoFilename = outputDir / ("component-%d-consensuses.txt" % count)
        if verbose:
            print(
                "      Saving component %d consensus FASTA to %s\n"
                "      Saving component %d consensus info to %s"
                % (count, consensusFilename, count, infoFilename)
            )
        with open(consensusFilename, "w") as consensusFp, open(
            infoFilename, "w"
        ) as infoFp:
            # First write the reference sequence for this component.
            (reference,) = list(
                FastaReads(outputDir / ("reference-component-%d.fasta" % count))
            )
            print(reference.toString("fasta"), file=consensusFp, end="")
            for i, cc in enumerate(self.consistentComponents, start=1):
                cc.saveConsensus(
                    i, self.offsets, referenceSequence, consensusFp, infoFp
                )

        # Write out an HTML table showing the identity between the various
        # component consensuses.
        identityTableFilename = outputDir / (
            "component-%d-consensuses-identity.html" % count
        )
        if verbose:
            print(
                "      Saving component %d consensus identity table to %s"
                % (count, identityTableFilename)
            )

        fastaIdentityTable(consensusFilename, identityTableFilename, verbose)


class ClusterAnalysis(ReadAnalysis):
    """
    Perform a clustered read alignment analysis for multiple infection
    detection.

    @param sampleName: The C{str} name of the sample whose reads are being
        analysed.
    @param alignmentFiles: A C{list} of C{str} names of SAM/BAM alignment
        files. These files should have mappings of reads to all references
        (as produced using the --all option when running bowtie2).
    @param referenceGenomeFiles: A C{list} of C{str} names of FASTA files
        containing reference genomes.
    @param outputDir: The C{Path} to the output directory.
    @param referenceIds: The C{str} sequence ids whose alignment should be
        analyzed. All ids must be present in the C{referenceGenomes} files.
        One of the SAM/BAM files given using C{alignmentFiles} should have an
        alignment against the given argument. If omitted, all references that
        are aligned to in the given BAM/SAM files will be analyzed.
    @param maxClusterDist: The C{float} inter-cluster distance above which
        clusters are too different to merge.
    @param minReads: The C{int} minimum number of reads that must cover a
        location for it to be considered significant.
    @param homogeneousCutoff: If the most common nucleotide at a location
        occurs more than this C{float} fraction of the time (i.e., amongst all
        reads that cover the location) then the locaion will be considered
        homogeneous and therefore uninteresting.
    @param plotSAM: If C{True} save plots of where reads lie on each reference
        genome (can be slow).
    @param alternateNucleotideMinFreq: The C{float} frequency that an
        alternative nucleotide (i.e., not the one chosen for the consensus)
        must have in order to be selected for the alternate consensus.
    @param minCCIdentity: The C{float} minimum nucleotide identity that a consistent
        component must have with a reference in order to contribute to the
        consensus being made against the reference.
    @param noCoverageStrategy: The approach to use when making a consensus if there are
        no reads covering a site. A value of 'N' means to use an ambigous N nucleotide
        code, whereas a value fo 'reference' means to take the base from the reference
        sequence.

    @param saveReducedFASTA: If C{True}, write out a FASTA file of the original
        input but with just the significant locations.
    @param verbose: The C{int} verbosity level. Use C{0} for no output.
    """

    DEFAULT_MAX_CLUSTER_DIST = 0.2
    ALTERNATE_NUCLEOTIDE_MIN_FREQ_DEF = 0.15
    MIN_CC_IDENTITY_DEFAULT = 0.7

    def __init__(
        self,
        readAnalysis,
        maxClusterDist=DEFAULT_MAX_CLUSTER_DIST,
        alternateNucleotideMinFreq=ALTERNATE_NUCLEOTIDE_MIN_FREQ_DEF,
        minCCIdentity=MIN_CC_IDENTITY_DEFAULT,
        noCoverageStrategy="N",
        plotSAM=False,
    ):
        self.readAnalysis = readAnalysis
        self.maxClusterDist = maxClusterDist
        self.alternateNucleotideMinFreq = alternateNucleotideMinFreq
        self.minCCIdentity = minCCIdentity
        assert noCoverageStrategy in {"N", "reference"}
        self.noCoverageStrategy = noCoverageStrategy

    def analyzeReference(self, reference):
        """
        Analyze a reference.

        @param reference: A C{Reference} instance to analyze.
        """

        components = self.findConnectedComponents(reference)
        self.saveComponentFasta(reference, components)
        self.summarize(reference, components)
        self.saveReferenceComponents(reference, components)
        self.saveComponentConsensuses(reference, components)

        (
            consensusRead,
            unwantedReads,
            wantedCcReadCount,
            consensusReadCountAtOffset,
            consensusWantedReadsBaseCountAtOffset,
        ) = self.saveClosestReferenceConsensus(
            reference, components, reference.paddedSAM.referenceInsertions
        )

        reference.consensusRead = consensusRead

        self.saveAlternateConsensus(reference, consensusRead)
        self.plotConsistentComponents(reference, components)
        self.saveConsensusBaseFrequencyPlot(
            reference,
            consensusWantedReadsBaseCountAtOffset,
            consensusReadCountAtOffset,
        )

    def findConnectedComponents(self, reference):
        """
        Find all connected components.

        @param reference: A C{Reference} instance to analyze.
        @return: A C{list} of C{Component} instances,
            sorted by component (the smallest offset is used for sorting
            so this gives the components from left to right along the
            reference genome.
        """
        significantReads = set(
            read for read in reference.alignedReads if read.significantOffsets
        )
        components = []
        filename = reference.outputDir / "cluster-analysis.txt"
        self.readAnalysis.report("    Saving clustering steps to", filename)

        with open(filename, "w") as fp:
            for count, component in enumerate(
                connectedComponentsByOffset(
                    significantReads, self.maxClusterDist, fp, self.readAnalysis.verbose
                ),
                start=1,
            ):
                components.append(component)

        # Sanity check: The significantReads set should be be empty
        # following the above processing by connectedComponentsByOffset.
        assert len(significantReads) == 0
        return sorted(components)

    def saveComponentFasta(self, reference, components):
        """
        Save FASTA for each component.

        @param outputDir: The C{Path} to the output directory.
        """
        self.readAnalysis.report("    Saving component FASTA")
        for count, component in enumerate(components, start=1):
            component.saveFasta(reference.outputDir, count, self.readAnalysis.verbose)

    def saveReferenceComponents(self, reference, components):
        """
        Save a FASTA file for the reference containing just the offsets for
        all connected components.

        @param reference: A C{Reference} instance.
        @param components: A C{list} of C{Component} instances.
        """
        for count, component in enumerate(components, start=1):
            filename = reference.outputDir / ("reference-component-%d.fasta" % count)
            self.readAnalysis.report(
                "    Saving reference component %d to %s" % (count, filename)
            )
            read = Read(reference.id + "-component-%d" % count, reference.read.sequence)

            Reads([read]).filter(keepSites=component.offsets).save(filename)

    def plotConsistentComponents(self, reference, components):
        """
        Make a plot of all consistent connected components.

        @param reference: A C{Reference} instance.
        @param components: A C{list} of C{Component} instances.
        """
        filename = reference.outputDir / "consistent-components-plot.html"
        self.readAnalysis.report(
            "    Plotting consistent connected components to", filename
        )
        infoFilename = reference.outputDir / "consistent-components-plot.txt"
        self.readAnalysis.report(
            "    Writing consistent connected component info to", infoFilename
        )

        plotConsistentComponents(
            reference,
            components,
            filename,
            infoFilename,
            titleFontSize=17,
            axisFontSize=15,
            title=(
                f"Consistent connected components when mapping "
                f"{len(reference.alignedReads)} reads from "
                f"{self.readAnalysis.sampleName} against "
                f"{reference.id}<br>from file "
                f"{self.readAnalysis.shortAlignmentFilename[reference.alignmentFile]}"
            ),
        )

    def saveClosestReferenceConsensus(self, reference, components, referenceInsertions):
        """
        Calculate and save the best consensus to a reference genome.

        @param reference: A C{Reference} instance.
        @param components: A C{list} of C{Component} instances.
        @param referenceInsertions: A C{dict} keyed by read id (the read
            that would cause a reference insertion). The values are lists
            of 2-tuples, with each 2-tuple containing an offset into the
            reference sequence and the C{str} of nucleotide that would be
            inserted starting at that offset.
        @return: A tuple of (consensus, unwantedReads, wantedCcReadCount,
                 wantedReadsCountAtOffset, wantedReadsBaseCountAtOffset).
        """

        def ccMatchCount(cc, reference, drawFp, drawMessage):
            """
            Count the matches between a consistent component and a reference
            genome.

            @param cc: A C{ConsistentComponent} instance.
            @param reference: A C{Reference} instance.
            @param drawFp: A file pointer to write information about draws (ifany) to.
            @param drawMessage: A C{str} message to write to C{drawFp}. If the
                string contains '%(baseCounts)s' that will be replaced by a
                string representation of the base counts (in C{counts})
                obtained from C{baseCountsToStr}. If not, the base count info
                will be printed after the message.
            @return: The C{int} count of bases that match the reference for the offsets
                covered by the consistent component.
            """
            referenceSequence = reference.read.sequence
            nucleotides = cc.nucleotides
            count = 0
            for offset in nucleotides:
                message = (
                    drawMessage + f" location {offset + 1}: base counts %(baseCounts)s."
                )
                referenceBase = referenceSequence[offset]
                componentBase = commonest(
                    nucleotides[offset],
                    referenceBase,
                    drawFp=drawFp,
                    drawMessage=message,
                )
                count += componentBase == referenceBase
            return count

        def scoreCcs(reference, component, fp):
            """
            Score the consistent components in the given C{Component}
            instance according to how well they match the passed reference.

            @param reference: A C{Reference} instance.
            @param component: A C{Component} instance.
            @param fp: A file pointer to write information to.
            @return: A C{list} of C{tuples}, each containing the score (nucleotide
                identity fraction against the reference), the index of the
                component, and the connected component.
            """
            result = []
            for index, cc in enumerate(component.consistentComponents):
                matchCount = ccMatchCount(
                    cc,
                    reference,
                    fp,
                    f"    Consistent component {index + 1} base draw",
                )
                score = matchCount / len(cc.nucleotides)
                print(
                    f"  Consistent component {index + 1} ({len(cc.reads)} reads) has "
                    f"{matchCount} exact matches with the reference, out of the "
                    f"{len(cc.nucleotides)} significant offsets it covers "
                    f"({score * 100.0:.2f}).",
                    file=fp,
                )
                result.append((score, index, cc))

            return result

        def partitionCcs(scoredCcs):
            """
            Partition the consistent components into high- (nucleotide identity
            with the reference) and low-scoring subsets.

            @return: A 2-C{tuple} containing two sets. Each set has either the
                high-scoring (i.e., high nucleotide identity fraction) or low-scoring
                consistent components, each with their index.
            """
            highScoringCCs = set()
            lowScoringCCs = set()

            for score, ccIndex, cc in scoredCcs:
                if score < self.minCCIdentity:
                    prefix = "Not I"
                    lowScoringCCs.add((ccIndex, cc))
                else:
                    prefix = "I"
                    highScoringCCs.add((ccIndex, cc))
                print(
                    "  %sncorporating nucleotides from consistent "
                    "component %d (%d reads, score %.2f, covering %d "
                    "locations (%d still undecided in consensus)) to "
                    "consensus."
                    % (
                        prefix,
                        ccIndex + 1,
                        len(cc.reads),
                        score,
                        len(cc.nucleotides),
                        len(set(cc.nucleotides) - offsetsDone),
                    ),
                    file=infoFp,
                )

            return highScoringCCs, lowScoringCCs

        referenceSequence = reference.read.sequence
        referenceLength = len(reference.read)
        fields = reference.read.id.split(maxsplit=1)
        referenceIdRest = "" if len(fields) == 1 else " " + fields[1]

        infoFile = reference.outputDir / "reference-consensus.txt"
        self.readAnalysis.report(
            "    Saving closest consensus to reference info to", infoFile
        )

        with open(infoFile, "w") as infoFp:
            print("Building consensus at significant offsets.", file=infoFp)
            consensus = [None] * referenceLength
            offsetsDone = set()
            wantedReads = set()
            unwantedReads = set()

            for componentCount, component in enumerate(components, start=1):
                print(
                    f"\nExamining component {componentCount} with "
                    f"{len(component.offsets)} locations: %s"
                    + commas(offset + 1 for offset in component.offsets),
                    file=infoFp,
                )
                scoredCcs = scoreCcs(reference, component, infoFp)

                highScoringCCs, lowScoringCCs = partitionCcs(scoredCcs)

                consistentComponent = ConsistentComponent()
                for _, cc in highScoringCCs:
                    consistentComponent.update(cc.reads)
                    wantedReads.update(cc.reads)

                for _, cc in lowScoringCCs:
                    unwantedReads.update(cc.reads)

                # Use the high-scoring CCs to determine consensus bases for the
                # significant offsets they collectively cover.

                for offset in sorted(consistentComponent.offsets):
                    assert offset not in offsetsDone
                    offsetsDone.add(offset)
                    assert consensus[offset] is None
                    consensus[offset] = consistentComponent.consensusBase(
                        offset, referenceSequence, infoFp
                    )
                    consensusBase = consensus[offset]
                    referenceBase = referenceSequence[offset]

                    mismatch = (
                        ""
                        if consensusBase == referenceBase
                        else (
                            f" (mismatch: reference has {referenceBase}, selected "
                            f"consisted component reads consensus has {consensusBase})"
                        )
                    )

                    print(
                        "    Location %d: %s from nucleotides %s%s"
                        % (
                            offset + 1,
                            consensusBase,
                            consistentComponent.nucleotides[offset].baseCountsToStr(),
                            mismatch,
                        ),
                        file=infoFp,
                    )

            # All components are now processed.

            # Get the base counts at each offset, from the full set of reads minus those
            # we're not using.
            (wantedReadsCountAtOffset, wantedReadsBaseCountAtOffset, _) = analyzeOffets(
                referenceLength, set(reference.alignedReads) - unwantedReads
            )
            remainingOffsets = sorted(set(range(referenceLength)) - offsetsDone)

            print(
                "\nAttempting to add bases from %d non-significant "
                "consensus locations, using all reads, EXCEPT those "
                "belonging to unused consistent components:" % len(remainingOffsets),
                file=infoFp,
            )

            uncoveredOffsetCount = 0
            for offset in remainingOffsets:
                baseCount = wantedReadsBaseCountAtOffset[offset]
                if baseCount:
                    referenceBase = reference.read.sequence[offset]
                    base = commonest(
                        baseCount,
                        referenceBase,
                        drawFp=infoFp,
                        drawMessage=(
                            "    WARNING: consensus base count draw at location %d"
                            % (offset + 1)
                        )
                        + " %(baseCounts)s.",
                    )
                    print(
                        "  Location %d: %s from nucleotides %s"
                        % (offset + 1, base, baseCountsToStr(baseCount)),
                        file=infoFp,
                        end="",
                    )

                    if base == referenceBase:
                        print(file=infoFp)
                    else:
                        print(
                            " (mismatch: reference has %s)" % referenceBase, file=infoFp
                        )
                else:
                    # An offset with no covering reads. Depending on the value of
                    # self.noCoverageStrategy, we either set this to "N" to indicate
                    # that there is (probably) a base there but that we don't have
                    # well-matching reads that cover the region, or we set it to the
                    # reference base.
                    uncoveredOffsetCount += 1
                    base = (
                        "N"
                        if self.noCoverageStrategy == "N"
                        else referenceSequence[offset]
                    )

                    print(
                        "  Location %d: %s NOT COVERED IN CONSENSUS%s"
                        % (
                            offset + 1,
                            base,
                            "."
                            if self.noCoverageStrategy == "N"
                            else " (set from reference).",
                        ),
                        file=infoFp,
                    )

                assert consensus[offset] is None
                consensus[offset] = base
                assert offset not in offsetsDone
                offsetsDone.add(offset)

            coveredOffsetCountOrig = sum(
                bool(counts) for counts in reference.baseCountAtOffset
            )
            coveredOffsetCountConsensus = referenceLength - uncoveredOffsetCount

            print(
                "\n"
                f"  Coverage based on all reads in original BAM file: "
                f"{pct(coveredOffsetCountOrig, referenceLength)}.",
                file=infoFp,
            )
            print(
                f"  Coverage in consensus "
                f"{pct(coveredOffsetCountConsensus, referenceLength)}.",
                file=infoFp,
            )
            if coveredOffsetCountConsensus != coveredOffsetCountOrig:
                assert coveredOffsetCountConsensus < coveredOffsetCountOrig
                coveredOffsetCountDiff = (
                    coveredOffsetCountOrig - coveredOffsetCountConsensus
                )
                print(
                    f"  The consensus has "
                    f"{pct(coveredOffsetCountDiff, referenceLength)} "
                    f"fewer base{s(coveredOffsetCountDiff)} covered.",
                    file=infoFp,
                )

            # Sanity check: make sure we processed all offsets.
            assert offsetsDone == set(range(referenceLength))

            consensusId = "%s-consensus%s" % (reference.id, referenceIdRest)

            consensus = Read(consensusId, "".join(consensus))

            # Print details of the match of the consensus to the reference.
            match = compareDNAReads(reference.read, consensus, matchAmbiguous=False)
            print("\nOVERALL match with reference:", file=infoFp)
            print("  Strict:", file=infoFp)
            print(
                matchToString(match, reference.read, consensus, indent="    "),
                file=infoFp,
            )

            match = compareDNAReads(reference.read, consensus, matchAmbiguous=True)
            print("  Non-strict:", file=infoFp)
            print(
                matchToString(match, reference.read, consensus, indent="    "),
                file=infoFp,
            )

            # Print any insertions to the reference.
            wantedReadsWithInsertions = set(referenceInsertions) & (
                set(reference.alignedReads) - unwantedReads
            )
            if wantedReadsWithInsertions:
                print(
                    "\nReference insertions present in %d read%s:"
                    % (
                        len(wantedReadsWithInsertions),
                        s(len(wantedReadsWithInsertions)),
                    ),
                    file=infoFp,
                )
                nucleotides = defaultdict(Counter)
                for readId in wantedReadsWithInsertions:
                    for offset, sequence in referenceInsertions[readId]:
                        for index, base in enumerate(sequence):
                            nucleotides[offset + index][base] += 1
                print(nucleotidesToStr(nucleotides, prefix="  "), file=infoFp)
            else:
                print("\nReference insertions: none.", file=infoFp)

        filename = reference.outputDir / "reference-consensus.fasta"
        self.readAnalysis.report("    Saving consensus to", filename)
        Reads([consensus]).save(filename)

        wantedCcReadCount = 0
        filename = reference.outputDir / "cc-wanted.fastq"
        with open(filename, "w") as fp:
            for wantedCcRead in wantedReads:
                alignment = wantedCcRead.alignment
                if not (alignment.is_secondary or alignment.is_supplementary):
                    wantedCcReadCount += 1
                    print(
                        Read(
                            alignment.query_name,
                            alignment.query_sequence,
                            alignmentQuality(alignment),
                        ).toString("fastq"),
                        end="",
                        file=fp,
                    )
        self.readAnalysis.report(
            f"    Saved {wantedCcReadCount} read{s(wantedCcReadCount)} wanted in "
            f"consistent components to {str(filename)!r}."
        )

        return (
            consensus,
            unwantedReads,
            wantedCcReadCount,
            wantedReadsCountAtOffset,
            wantedReadsBaseCountAtOffset,
        )

    def saveAlternateConsensus(self, reference, consensusRead):
        """
        Calculate and save an alternate consensus to a reference genome.

        @param reference: A C{Reference} instance.
        @consensusRead: The C{dark.reads.Read} consensus sequence to calculate
            an alternative to.
        @return: An alternate consensus C{dark.reads.Read} instance.
        """
        filename = reference.outputDir / "reference-alternate-consensus.txt"
        self.readAnalysis.report(
            f"    Writing alternate consensus info to {str(filename)!r}"
        )
        alternateConsensus = []

        fields = reference.read.id.split(maxsplit=1)
        referenceIdRest = "" if len(fields) == 1 else " " + fields[1]

        with open(filename, "w") as infoFp:
            print(
                "The three nucleotides shown on each line are for the "
                "reference, the consensus (made by clustering) and the "
                "alternate nucleotide.\nThese are followed by up to two "
                "asterisks: the first if the alternate nucleotide does not "
                "agree with the reference, the second if it does not agree "
                "with the consensus.",
                file=infoFp,
            )
            for offset in range(len(reference.read)):
                referenceBase = reference.read.sequence[offset]
                consensusBase = consensusRead.sequence[offset]
                baseCount = reference.baseCountAtOffset[offset]
                if baseCount:
                    if len(baseCount) == 1:
                        # Only one nucleotide was found at this location.  The reference
                        # doesn't necessarily agree with the consensus here, since the
                        # aligned reads may have had a different base at this site.
                        base = consensusBase
                    else:
                        # Find the nucleotide with the highest count that is not the
                        # consensus sequence nucleotide.
                        orderedCounts = baseCount.most_common()
                        alternateBase, alternateCount = [
                            x for x in orderedCounts if x[0] != consensusBase
                        ][0]

                        # Check that we found a base that's not the consensus base.
                        assert alternateBase != consensusBase

                        # If the frequency of the alternate base is high enough, go with
                        # it. Else take the base from the original consensus.
                        alternateFraction = (
                            alternateCount / reference.readCountAtOffset[offset]
                        )
                        if (
                            alternateCount > 1
                            and alternateFraction > self.alternateNucleotideMinFreq
                        ):
                            base = alternateBase
                        else:
                            base = consensusBase

                    agreeWithReference = referenceBase == base
                    agreeWithConsensus = consensusBase == base
                    print(
                        "Location %d: %s %s %s %s %s nucleotides %s"
                        % (
                            offset + 1,
                            referenceBase,
                            consensusBase,
                            base,
                            " " if agreeWithReference else "*",
                            " " if agreeWithConsensus else "*",
                            baseCountsToStr(baseCount),
                        ),
                        file=infoFp,
                    )
                else:
                    # The reads did not cover this offset.
                    base = "-"
                    print("  Location %d: -" % (offset + 1), file=infoFp)

                alternateConsensus.append(base)

            alternateConsensusId = "%s-alternate-consensus%s" % (
                reference.id,
                referenceIdRest,
            )

            alternateConsensusRead = Read(
                alternateConsensusId, "".join(alternateConsensus)
            )

            # Print details of the match of the alternate consensus to the
            # reference.
            match = compareDNAReads(reference.read, alternateConsensusRead)
            print("\nAlternate consensus match with reference:", file=infoFp)
            print(
                matchToString(
                    match, reference.read, alternateConsensusRead, indent="  "
                ),
                file=infoFp,
            )

            # Print details of the match of the alternate consensus to the
            # original consensus.
            match = compareDNAReads(consensusRead, alternateConsensusRead)
            print("\nAlternate consensus match with original consensus:", file=infoFp)
            print(
                matchToString(
                    match, consensusRead, alternateConsensusRead, indent="  "
                ),
                file=infoFp,
            )

            # Print details of the match of the original consensus to the
            # reference.
            match = compareDNAReads(reference.read, consensusRead)
            print("\nOriginal consensus match with reference:", file=infoFp)
            print(
                matchToString(match, reference.read, consensusRead, indent="  "),
                file=infoFp,
            )

        filename = reference.outputDir / "reference-alternate-consensus.fasta"
        self.readAnalysis.report("    Saving alternate consensus FASTA to", filename)
        Reads([alternateConsensusRead]).save(filename)

        return alternateConsensusRead

    def saveConsensusBaseFrequencyPlot(
        self,
        reference,
        baseCountAtOffset,
        readCountAtOffset,
    ):
        """
        Make a plot of the sorted base frequencies for the consensus.

        @param reference: A C{Reference} instance.
        @param baseCountAtOffset: A C{list} of C{Counter} instances giving
            the count of each nucleotide at each genome offset.
        @param readCountAtOffset: A C{list} of C{int} counts of the total
            number of reads at each genome offset (i.e., just the sum of the
            values in C{baseCountAtOffset})
        """
        filename = reference.outputDir / "consensus-base-frequencies.html"
        self.readAnalysis.report(
            "    Writing consensus base frequency plot to", filename
        )

        title = (
            f"Significant sites base frequencies when mapping "
            f"{len(reference.alignedReads)} reads from {self.readAnalysis.sampleName} "
            f"against the {reference.id} consensus (length {len(reference.read)})."
        )

        significantOffsets = list(
            findSignificantOffsets(
                baseCountAtOffset,
                readCountAtOffset,
                self.readAnalysis.minReads,
                self.readAnalysis.homogeneousCutoff,
            )
        )

        plotBaseFrequencies(
            len(reference.read),
            significantOffsets,
            baseCountAtOffset,
            readCountAtOffset,
            filename,
            title=title,
            minReads=self.readAnalysis.minReads,
            homogeneousCutoff=self.readAnalysis.homogeneousCutoff,
            histogram=False,
            show=False,
        )

    def saveComponentConsensuses(self, reference, components):
        """
        Write out a component consensus sequence.

        @param reference: A C{Reference} instance.
        @param components: A C{list} of C{Component} instances.
        """
        self.readAnalysis.report("    Saving component consensuses")
        for count, component in enumerate(components, start=1):
            component.saveConsensuses(
                reference.outputDir,
                count,
                reference.read.sequence,
                self.readAnalysis.verbose,
            )

    def summarize(self, reference, components):
        """
        Write out an analysis summary.

        @param reference: A C{Reference} instance.
        @param components: A C{list} of C{Component} instances.
        """
        filename = reference.outputDir / "component-summary.txt"
        self.readAnalysis.report(f"    Writing analysis summary to {quoted(filename)}.")

        with open(filename, "w") as fp:
            print(
                "Read %d aligned reads of length %d. Found %d significant locations."
                % (
                    len(reference.alignedReads),
                    len(reference.read),
                    len(reference.significantOffsets),
                ),
                file=fp,
            )

            print(
                "Reads were assigned to %d connected components:" % len(components),
                file=fp,
            )

            totalReads = 0
            for count, component in enumerate(components, start=1):
                filename = reference.outputDir / ("component-%d.txt" % count)
                self.readAnalysis.report(
                    f"    Writing component {count} summary to {quoted(filename)}."
                )
                with open(filename, "w") as fp2:
                    component.summarize(fp2, count, reference.read.sequence)

                componentCount = len(component)
                offsets = component.offsets
                totalReads += componentCount
                print(
                    "\nConnected component %d: %d reads, covering %d offsets "
                    "(%d to %d)"
                    % (count, componentCount, len(offsets), min(offsets), max(offsets)),
                    file=fp,
                )

                ccCounts = sorted(
                    map(
                        len,
                        (cc.reads for cc in component.consistentComponents),
                    ),
                    reverse=True,
                )
                if len(ccCounts) > 1:
                    print(
                        "  largest two consistent component size ratio %.2f"
                        % (ccCounts[0] / ccCounts[1]),
                        file=fp,
                    )

                for j, cc in enumerate(component.consistentComponents, start=1):
                    print(
                        "  consistent sub-component %d: read count %d, covered offset "
                        "count %d."
                        % (
                            j,
                            len(cc.reads),
                            len(cc.nucleotides),
                        ),
                        file=fp,
                    )

            print(
                "\nIn total, %d reads were assigned to components." % totalReads,
                file=fp,
            )

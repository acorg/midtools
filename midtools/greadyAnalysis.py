from os.path import join
from pysam import AlignmentFile

from dark.reads import Read, Reads
from dark.sam import samfile

from midtools.analysis import ReadAnalysis
from midtools.offsets import OffsetBases
from midtools.pqueue import PriorityQueue
from midtools.utils import s


class GreadyAnalysis(ReadAnalysis):
    """
    Perform a greedy read alignment analysis for multiple infection detection.

    @param alignmentFiles: A C{list} of C{str} names of SAM/BAM alignment
        files.
    @param referenceGenomeFiles: A C{list} of C{str} names of FASTA files
        containing reference genomes.
    @param referenceIds: The C{str} sequence ids whose alignment should be
        analyzed. All ids must be present in the C{referenceGenomes} files.
        One of the SAM/BAM files given using C{alignmentFiles} should have an
        alignment against the given argument. If omitted, all references that
        are aligned to in the given BAM/SAM files will be analyzed.
    @param cutoff: Reads with a significant offset C{float} mismatch fraction
        greater than this will not be put into the consensus alignment.
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

    DEFAULT_AGREEMENT_THRESHOLD = 0.5

    def __init__(
        self,
        alignmentFiles,
        referenceGenomeFiles,
        referenceIds=None,
        cutoff=0.5,
        outputDir=None,
        minReads=ReadAnalysis.DEFAULT_MIN_READS,
        homogeneousCutoff=ReadAnalysis.DEFAULT_HOMOGENEOUS_CUTOFF,
        plotSAM=False,
        saveReducedFASTA=False,
        verbose=0,
    ):
        super().__init__(
            alignmentFiles,
            referenceGenomeFiles,
            referenceIds=referenceIds,
            outputDir=outputDir,
            minReads=minReads,
            homogeneousCutoff=homogeneousCutoff,
            plotSAM=plotSAM,
            saveReducedFASTA=saveReducedFASTA,
            verbose=verbose,
        )

        self.cutoff = cutoff

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
            referenceId, alignmentFile, outputDir
        )

        if analysis:
            (
                genomeLength,
                alignedReads,
                readCountAtOffset,
                baseCountAtOffset,
                readsAtOffset,
                significantOffsets,
                samFilter,
                paddedSAM,
            ) = analysis
        else:
            return

        insignificantOffsets = set(range(genomeLength)) - set(significantOffsets)

        reference = self.referenceGenomes[referenceId]
        referenceSequence = reference.sequence

        consensus = []
        for base in referenceSequence:
            ob = OffsetBases()
            ob.incorporateBase(base)
            consensus.append(ob)

        readQueue = PriorityQueue()
        self.updatePriorityQueue(readQueue, alignedReads, consensus, significantOffsets)

        consensusFilename = join(outputDir, "reference-consensus.sam")
        nonConsensusFilename = join(outputDir, "reference-non-consensus.sam")
        self.report("    Writing consensus SAM to", consensusFilename)
        self.report("    Writing non-consensus SAM to", nonConsensusFilename)

        with samfile(alignmentFile) as sam:
            consensusAlignment = AlignmentFile(
                consensusFilename, mode="w", template=sam
            )
            nonConsensusAlignment = AlignmentFile(
                nonConsensusFilename, mode="w", template=sam
            )

        # Reads with no significant offsets get written to both output files.
        readsWithNoSignificantOffsetsCount = 0
        for read in alignedReads:
            if not read.significantOffsets:
                readsWithNoSignificantOffsetsCount += 1
                consensusAlignment.write(read.alignment)
                nonConsensusAlignment.write(read.alignment)

                for offset in insignificantOffsets:
                    base = read.base(offset)
                    if base is not None:
                        consensus[offset].incorporateBase(base)

        self.report(
            "    %d read%s did not overlap any significant offsets"
            % (
                readsWithNoSignificantOffsetsCount,
                s(readsWithNoSignificantOffsetsCount),
            )
        )

        readsMatchingConsensusCount = readsNotMatchingConsensusCount = 0
        cutoff = self.cutoff
        while readQueue:
            mismatchFraction, _ = readQueue.lowestPriority()
            read = readQueue.pop()
            if mismatchFraction <= cutoff:
                # We want this read. Incorporate it into the consensus.
                readsMatchingConsensusCount += 1
                consensusAlignment.write(read.alignment)
                affectedReads = set()
                for offset in read.significantOffsets:
                    readBase = read.base(offset)
                    consensus[offset].incorporateBase(readBase)
                    for readAtOffset in readsAtOffset[offset]:
                        if readAtOffset in readQueue:
                            affectedReads.add(readAtOffset)
                self.updatePriorityQueue(
                    readQueue, affectedReads, consensus, significantOffsets
                )
            else:
                readsNotMatchingConsensusCount += 1
                nonConsensusAlignment.write(read.alignment)

        consensusAlignment.close()
        nonConsensusAlignment.close()

        self.report(
            "    %d read%s matched the consensus, %d did not."
            % (
                readsMatchingConsensusCount,
                s(readsMatchingConsensusCount),
                readsNotMatchingConsensusCount,
            )
        )

        # Remove the reference bases from the consensus.
        for offset, base in enumerate(referenceSequence):
            consensus[offset].unincorporateBase(base)

        consensusInfoFilename = join(outputDir, "reference-consensus.txt")
        self.report("    Writing consensus info to", consensusInfoFilename)

        with open(consensusInfoFilename, "w") as fp:
            consensusSequence = []
            for offset in range(genomeLength):
                # Take a copy of the commonest set because we may pop from
                # it below.
                commonest = set(consensus[offset].commonest)
                referenceBase = referenceSequence[offset]

                if len(commonest) > 1:
                    nucleotides = " Nucleotides: %s" % (
                        consensus[offset].baseCountsToStr()
                    )
                else:
                    nucleotides = ""

                if referenceBase in commonest:
                    consensusBase = referenceBase
                else:
                    if len(commonest) == 1:
                        # Nothing in the included reads covers this offset.
                        consensusBase = "-"
                    elif len(commonest) > 1:
                        # Report a draw (in which the reference base is not
                        # included and so cannot be used to break the draw).
                        commonest.pop()
                    else:
                        consensusBase = commonest.pop()

                consensusSequence.append(consensusBase)

                mismatch = (
                    ""
                    if referenceBase == consensusBase
                    else (" Mismatch (reference has %s)" % referenceBase)
                )

                print(
                    "%d: %s%s%s" % (offset + 1, consensusBase, mismatch, nucleotides),
                    file=fp,
                )

        consensusRead = Read(
            "gready-consensus-%s" % referenceId, "".join(consensusSequence)
        )
        consensusFilename = join(outputDir, "reference-consensus.fasta")
        self.report("    Writing gready consensus info to", consensusFilename)
        Reads([consensusRead]).save(consensusFilename)

        return {
            "consensusRead": consensusRead,
            "significantOffsets": significantOffsets,
        }

    def updatePriorityQueue(self, readQueue, reads, consensus, significantOffsets):
        """
        Put all aligned reads that cover significant offsets into the priority
        queue, giving each a score corresponding to how well it matches the
        consensus.

        @param readQueue: A C{pqueue.PriorityQueue} instance.
        @param reads: A list of C{AlignedRead} instances.
        @param consensus: A C{list} of C{OffsetBases} instances holding the
            base counts for all offsets.
        """
        significantOffsetCount = len(significantOffsets)

        for read in reads:
            if read.significantOffsets:
                matches = 0
                for offset, base in read.significantOffsets.items():
                    matches += base in consensus[offset].commonest

                # Lower priority is better.
                mismatchFraction = 1.0 - (matches / len(read.significantOffsets))
                significantOffsetsMissing = significantOffsetCount - len(
                    read.significantOffsets
                )

                priority = (mismatchFraction, significantOffsetsMissing)

                readQueue.add(read, priority)

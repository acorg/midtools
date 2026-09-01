#!/usr/bin/env python

from itertools import chain

from midtools.analysis import ReadAnalysis
from midtools.clusterAnalysis import ClusterAnalysis
from midtools.options import addAnalysisCommandLineOptions

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Make consensus sequences by aligning reads to references "
            "and finding which reads agree and disagree with one another "
            "using clustering."
        ),
    )

    addAnalysisCommandLineOptions(parser)

    parser.add_argument(
        "--maxClusterDist",
        type=float,
        default=ClusterAnalysis.DEFAULT_MAX_CLUSTER_DIST,
        help=(
            "Clustering will be stopped once the minimum distance between "
            "remaining clusters exceeds this value."
        ),
    )

    parser.add_argument(
        "--alternateNucleotideMinFreq",
        type=float,
        default=ClusterAnalysis.ALTERNATE_NUCLEOTIDE_MIN_FREQ_DEF,
        help=(
            "The (0.0 to 1.0) frequency that an alternative nucleotide "
            "(i.e., not the one chosen for a consensus) must have in order "
            "to be selected for the alternate consensus."
        ),
    )

    parser.add_argument(
        "--minCCIdentity",
        type=float,
        default=ClusterAnalysis.MIN_CC_IDENTITY_DEFAULT,
        help=(
            "The minimum nucleotide identity fraction [0.0, 1.0] that a consistent "
            "component must have with a reference in order to contribute to the "
            "consensus being made against the reference."
        ),
    )

    parser.add_argument(
        "--noCoverageStrategy",
        default="N",
        choices=("N", "reference"),
        help=(
            "The approach to use when making a consensus if there are no reads "
            "covering a site. A value of 'N' means to use an ambigous N nucleotide "
            "code, whereas a value of 'reference' means to take the base from the "
            "reference sequence."
        ),
    )

    args = parser.parse_args()

    referenceIds = (
        list(chain.from_iterable(args.referenceId)) if args.referenceId else None
    )

    # The logic of the below is a little hard to follow. There are three steps:
    #
    # 1. Make a read analysis with a 'run' method (which we don't call yet).
    #
    # 2. Set up a Cluster Analysis, giving it the read analysis (so the running cluster
    #    analysis can get some parameters and produce reporting information when the
    #    read analysis is run.
    #
    # 3. Call the read analysis 'run' method (see 1 above), passing it the function from
    #    the cluster analysis that can analyze a reference.
    #
    # Things are done in this (seemingly?) convoluted way because I implemented several
    # methods for doing the main work, of which the cluster analysis is just one. They
    # all needed to work in the same way, and so I separated the organizational things
    # (like making output directories and collecting final / overall results) into the
    # common read analysis class, and the code that does the actual analysis for an
    # input alignment file / reference pair. If you look at the 'run' method of the
    # ReadAnalysis class you'll see it's basically just a loop over alignment file /
    # reference pairs, each with some setup and then calling the cluster analysis
    # function. Then some final gathering of overall results.

    analysis = ReadAnalysis(
        args.sampleName,
        list(chain.from_iterable(args.alignmentFile)),
        list(chain.from_iterable(args.referenceGenome)),
        args.outputDir,
        referenceIds=referenceIds,
        minReads=args.minReads,
        homogeneousCutoff=args.homogeneousCutoff,
        plotSAM=args.plotSAM,
        saveReducedFASTA=args.saveReducedFASTA,
        verbose=args.verbose,
    )

    clusterAnalysis = ClusterAnalysis(
        analysis,
        maxClusterDist=args.maxClusterDist,
        alternateNucleotideMinFreq=args.alternateNucleotideMinFreq,
        minCCIdentity=args.minCCIdentity,
        noCoverageStrategy=args.noCoverageStrategy,
    )

    analysis.run(clusterAnalysis.analyzeReference)

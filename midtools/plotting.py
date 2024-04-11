import sys
from random import uniform
import plotly
from plotly.subplots import make_subplots
import plotly.graph_objs as go
import plotly.express as px
from operator import itemgetter
from json import dump
from itertools import cycle
from collections import Counter, defaultdict
from textwrap import wrap

from dark.dna import compareDNAReads
from dark.fasta import FastaReads

from midtools.entropy import entropy2, MAX_ENTROPY
from midtools.match import matchToString
from midtools.utils import s, baseCountsToStr

# The following import will fail for those that don't have our hbv repo.
# For now, just let the ImportError happen.
from pyhbv.genotype import getGenotype, genotypeKey
from pyhbv.samples import UNKNOWN, sampleIdKey


def plotSAM(
    samFilter,
    outfile,
    title="Reads",
    titleFontSize=18,
    axisFontSize=16,
    show=False,
    jitter=0.0,
):
    """
    Plot the alignments found in a SAM file.
    """
    referenceLengths = samFilter.referenceLengths()

    if len(set(referenceLengths.values())) == 1:
        _, referenceLength = referenceLengths.popitem()
    else:
        raise ValueError(
            "SAM/BAM file reference sequences lengths (%s) are not "
            "all identical." % ", ".join(map(str, sorted(referenceLengths)))
        )

    data = []

    count = 0
    for count, alignment in enumerate(samFilter.alignments(), start=1):
        referenceStart = alignment.reference_start
        score = alignment.get_tag("AS") + (
            0.0 if jitter == 0.0 else uniform(-jitter, jitter)
        )
        id_ = alignment.query_name
        data.append(
            go.Scatter(
                x=(referenceStart, referenceStart + alignment.reference_length),
                y=(score, score),
                text=(id_, id_),
                hoverinfo="text",
                mode="lines",
                showlegend=False,
            )
        )

    xaxis = {
        "title": "Genome location",
        "range": (0, referenceLength),
        "titlefont": {
            "size": axisFontSize,
        },
    }

    yaxis = {
        "title": "score",
        "titlefont": {
            "size": axisFontSize,
        },
    }

    title = "%s<br>%d reads%s" % (
        title,
        count,
        "" if jitter == 0.0 else (" (jitter %.2f)" % jitter),
    )

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            "size": titleFontSize,
        },
        hovermode="closest",
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def plotAllReferencesSAM(
    samFilter,
    outfile,
    sampleName,
    alignmentFile,
    hbv=False,
    jitter=0.0,
    titleFontSize=14,
    axisFontSize=12,
    show=False,
):
    """
    Plot the alignments found in a SAM file, even if the references are different
    lengths.

    @param sampleName: The C{str} name of the sample whose reads are being
        analysed.
    """
    referenceLengths = samFilter.referenceLengths()
    maxReferenceLength = max(referenceLengths.values())

    if hbv:
        referenceScaleFactor = dict(
            (id_, maxReferenceLength / length)
            for id_, length in referenceLengths.items()
        )
        referenceGenotype = dict(
            (id_, getGenotype(id_) or UNKNOWN) for id_ in referenceLengths
        )
    else:
        referenceScaleFactor = dict.fromkeys(referenceLengths, 1.0)
        referenceGenotype = dict((id_, id_) for id_ in referenceLengths)

    genotypes = sorted(set(referenceGenotype.values()), key=genotypeKey)
    legendRank = dict((genotype, i) for i, genotype in enumerate(genotypes))
    nGenotypes = len(genotypes)
    genotypeReferences = defaultdict(set)

    data = []
    inLegend = set()
    genotypeCount = Counter()

    count = 0
    for count, alignment in enumerate(samFilter.alignments(), start=1):
        referenceId = alignment.reference_name
        scaleFactor = referenceScaleFactor[referenceId]
        start = alignment.reference_start
        end = start + alignment.reference_length
        score = alignment.get_tag("AS") + (
            0.0 if jitter == 0.0 else uniform(-jitter, jitter)
        )
        genotype = referenceGenotype[referenceId]
        genotypeCount[genotype] += 1
        genotypeReferences[genotype].add(referenceId)
        text = f"{referenceId} Match {start + 1}-{end}, Read: {alignment.query_name}"
        data.append(
            go.Scatter(
                x=((start + 1) * scaleFactor, end * scaleFactor),
                y=(score, score),
                text=(text, text),
                hoverinfo="text",
                mode="lines",
                name=genotype,
                legendgroup=genotype,
                legendrank=legendRank[genotype],
                showlegend=genotype not in inLegend,
            )
        )
        inLegend.add(genotype)

    xaxis = {
        "title": "HBV genome location (scaled)" if hbv else "Genome location",
        "range": (0, maxReferenceLength),
        "titlefont": {
            "size": axisFontSize,
        },
    }

    yaxis = {
        "title": "Bowtie2 alignment score",
        "titlefont": {
            "size": axisFontSize,
        },
    }

    if hbv:
        genotypeReferencesDesc = (
            " Genotypes: "
            + "; ".join(
                (
                    f"<b>{gt}</b>: "
                    + ", ".join(sorted(genotypeReferences[gt], key=sampleIdKey))
                )
                for gt in genotypes
                if genotypeCount[gt]
            )
            + "."
        )
        sampleGenotypeDesc = f"(genotype {getGenotype(sampleName)}) "
    else:
        genotypeReferencesDesc = sampleGenotypeDesc = ""

    title = "<br>".join(
        wrap(
            f"Best-matched genotypes for {count} reads for {sampleName} "
            f"{sampleGenotypeDesc}from {alignmentFile}.{genotypeReferencesDesc}",
            width=175,
        )
    )

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            "size": titleFontSize,
        },
        hovermode="closest",
    )
    fig = go.Figure(data=data, layout=layout)

    # See "Color Sequences in Plotly Express" at
    # https://plotly.com/python/discrete-color/ for color sequences.
    colors = px.colors.qualitative.D3
    if nGenotypes > len(colors):
        print(
            f"WARNING: You have more genotypes ({nGenotypes}) than unique "
            f"colors ({len(colors)}). Some colors will be repeated.",
            file=sys.stderr,
        )

    genotypeColor = {}
    iterColors = cycle(colors)
    for genotype in genotypes:
        if genotypeCount[genotype]:
            genotypeColor[genotype] = next(iterColors)

    # Put the genotype read count into the legend labels and add the colors.
    fig.for_each_trace(
        lambda t: t.update(
            name=f"{t.name} ({genotypeCount[t.name]})",
            marker_color=genotypeColor[t.name],
        )
    )

    fig.update_layout(legend_title_text="Genotype (total reads)")

    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def _plotSortedMaxBaseFrequencies(
    significantOffsets,
    baseCountAtOffset,
    readCountAtOffset,
    outfile,
    title,
    histogram,
    show,
    titleFontSize,
    axisFontSize,
):
    """
    Plot the sorted maximum base frequency for each of the significant
    offsets.
    """
    frequencyInfo = []

    for offset in significantOffsets:
        count = readCountAtOffset[offset]

        sortedFreqs = [
            x / count for x in sorted(baseCountAtOffset[offset].values(), reverse=True)
        ]

        text = "site %d<br>" % (offset + 1) + ", ".join(
            "%s: %d" % (k, v) for k, v in baseCountAtOffset[offset].items()
        )

        frequencyInfo.append((sortedFreqs[0], text))

    # We don't have to sort if we're making a histogram, but we're expected
    # to return a sorted values list, so we sort unconditionally.
    frequencyInfo.sort(key=itemgetter(0))
    values = [freq for freq, _ in frequencyInfo]

    if histogram:
        data = [
            go.Histogram(x=values, histnorm="probability"),
        ]

        xaxis = {
            "title": "Significant site maximum nucleotide frequency",
            "range": (-0.05, 1.05),
            "titlefont": {
                "size": axisFontSize,
            },
        }

        yaxis = {
            "title": "Probability mass",
            "range": (0.0, 1.0),
            "titlefont": {
                "size": axisFontSize,
            },
        }
    else:
        data = [
            go.Scatter(
                x=list(range(1, len(significantOffsets) + 1)),
                y=values,
                mode="markers",
                showlegend=False,
                text=[text for _, text in frequencyInfo],
            ),
        ]

        xmargin = max(1, int(len(significantOffsets) * 0.01))
        xaxis = {
            "title": "Rank",
            "range": (-xmargin, len(significantOffsets) + xmargin),
            "titlefont": {
                "size": axisFontSize,
            },
        }

        yaxis = {
            "range": (0.0, 1.05),
            "title": "Significant site maximum nucleotide frequency",
            "titlefont": {
                "size": axisFontSize,
            },
        }

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            "size": titleFontSize,
        },
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)
    return frequencyInfo


def _plotBaseFrequenciesEntropy(
    significantOffsets,
    baseCountAtOffset,
    readCountAtOffset,
    outfile,
    title,
    histogram,
    show,
    titleFontSize,
    axisFontSize,
):
    """
    Plot the sorted entropy of base frequencies for each of the significant
    offsets.
    """
    entropyInfo = []

    for offset in significantOffsets:
        text = "site %d<br>" % (offset + 1) + ", ".join(
            "%s: %d" % (k, v) for k, v in baseCountAtOffset[offset].items()
        )

        entropyInfo.append((entropy2(list(baseCountAtOffset[offset].elements())), text))

    assert all([ent <= MAX_ENTROPY for ent, _ in entropyInfo])

    # We don't have to sort if we're making a histogram, but we're expected
    # to return a sorted values list, so we sort unconditionally.
    entropyInfo.sort(key=itemgetter(0))
    values = [ent for ent, _ in entropyInfo]

    if histogram:
        data = [go.Histogram(x=values, histnorm="probability")]

        xaxis = {
            "title": ("Significant site nucleotide frequency entropy (bits)"),
            "range": (-0.05, MAX_ENTROPY),
            "titlefont": {
                "size": axisFontSize,
            },
        }

        yaxis = {
            "title": "Probability mass",
            "range": (0.0, 1.0),
            "titlefont": {
                "size": axisFontSize,
            },
        }
    else:
        data = [
            go.Scatter(
                x=list(range(1, len(significantOffsets) + 1)),
                y=values,
                mode="markers",
                showlegend=False,
                text=[text for _, text in entropyInfo],
            ),
        ]

        xmargin = max(1, int(len(significantOffsets) * 0.01))
        xaxis = {
            "range": (-xmargin, len(significantOffsets) + xmargin),
            "title": "Rank",
            "titlefont": {
                "size": axisFontSize,
            },
        }

        yaxis = {
            "range": (-0.05, MAX_ENTROPY),
            "title": ("Significant site nucleotide frequency entropy (bits)"),
            "titlefont": {
                "size": axisFontSize,
            },
        }

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            "size": titleFontSize,
        },
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)
    return entropyInfo


def _plotBaseFrequenciesAllOffsets(
    genomeLength,
    significantOffsets,
    baseCountAtOffset,
    readCountAtOffset,
    outfile,
    title,
    show,
    titleFontSize,
    axisFontSize,
    yRange,
):
    """
    Plot the (sorted) base frequencies for each of the significant offsets.
    """

    # This function is currently unused. It plots all offsets, but the data become
    # invisible when there are under 350 significant offsets, so it's pretty useless
    # (you have to zoom to see any data).
    x = list(range(genomeLength))
    text = []
    freqs = [], [], [], []

    for offset in range(genomeLength):
        if offset in significantOffsets:
            count = readCountAtOffset[offset]

            sortedFreqs = [
                x / count
                for x in sorted(baseCountAtOffset[offset].values(), reverse=True)
            ]
            while len(sortedFreqs) < 4:
                sortedFreqs.append(0.0)

            for i, frequency in enumerate(sortedFreqs):
                freqs[i].append(frequency)

            text.append(
                ("site %d<br>" % (offset + 1))
                + ", ".join(
                    "%s: %d" % (k, v) for k, v in baseCountAtOffset[offset].items()
                )
            )
        else:
            for i in 0, 1, 2, 3:
                freqs[i].append(0.0)
            text.append("")

    data = [
        go.Bar(x=x, y=freqs[0], showlegend=False, text=text),
        go.Bar(x=x, y=freqs[1], showlegend=False),
        go.Bar(x=x, y=freqs[2], showlegend=False),
        go.Bar(x=x, y=freqs[3], showlegend=False),
    ]
    layout = go.Layout(
        barmode="stack",
        title=title,
        titlefont={
            "size": titleFontSize,
        },
        xaxis={
            "title": "Significant site index",
            "titlefont": {
                "size": axisFontSize,
            },
        },
        yaxis={
            "title": "Nucleotide frequency",
            "range": yRange,
            "titlefont": {
                "size": axisFontSize,
            },
        },
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def _plotBaseFrequencies(
    significantOffsets,
    baseCountAtOffset,
    readCountAtOffset,
    outfile,
    title,
    show,
    titleFontSize,
    axisFontSize,
    yRange,
):
    """
    Plot the (sorted) base frequencies for each of the significant offsets.
    """
    x = list(range(len(significantOffsets)))
    text = []
    freqs = [], [], [], []

    for offset in significantOffsets:
        count = readCountAtOffset[offset]

        sortedFreqs = [
            x / count for x in sorted(baseCountAtOffset[offset].values(), reverse=True)
        ]
        while len(sortedFreqs) < 4:
            sortedFreqs.append(0.0)

        for i, frequency in enumerate(sortedFreqs):
            freqs[i].append(frequency)

        text.append(
            ("site %d<br>" % (offset + 1))
            + ", ".join("%s: %d" % (k, v) for k, v in baseCountAtOffset[offset].items())
        )

    data = [
        go.Bar(x=x, y=freqs[0], showlegend=False, text=text),
        go.Bar(x=x, y=freqs[1], showlegend=False),
        go.Bar(x=x, y=freqs[2], showlegend=False),
        go.Bar(x=x, y=freqs[3], showlegend=False),
    ]
    layout = go.Layout(
        barmode="stack",
        title=title,
        titlefont={
            "size": titleFontSize,
        },
        xaxis={
            "title": "Significant site index",
            "titlefont": {
                "size": axisFontSize,
            },
        },
        yaxis={
            "title": "Nucleotide frequency",
            "range": yRange,
            "titlefont": {
                "size": axisFontSize,
            },
        },
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def plotBaseFrequencies(
    genomeLength,
    significantOffsets,
    baseCountAtOffset,
    readCountAtOffset,
    outfile,
    title=None,
    sampleName=None,
    valuesFile=None,
    minReads=5,
    homogeneousCutoff=0.9,
    sortOn=None,
    histogram=False,
    show=False,
    titleFontSize=12,
    axisFontSize=12,
    yRange=(0.0, 1.0),
):
    """
    Plot sorted base frequencies at signifcant sites.

    @param sampleName: The C{str} name of the sample whose reads are being
        analysed.
    """

    subtitle = (
        "<br>%d significant sites. Min %d read%s per site. "
        "%.2f homogeneity cutoff."
        % (len(significantOffsets), minReads, s(minReads), homogeneousCutoff)
    )

    if sortOn is None:
        title = title or "Base frequencies (sorted)"
        _plotBaseFrequencies(
            significantOffsets,
            baseCountAtOffset,
            readCountAtOffset,
            outfile,
            title + subtitle,
            show,
            titleFontSize,
            axisFontSize,
            yRange,
        )
    elif sortOn == "max":
        title = title or "Maximum base frequency"
        result = _plotSortedMaxBaseFrequencies(
            significantOffsets,
            baseCountAtOffset,
            readCountAtOffset,
            outfile,
            title + subtitle,
            histogram,
            show,
            titleFontSize,
            axisFontSize,
        )
    else:
        assert sortOn == "entropy", f"Unknown --sortOn value: {sortOn!r}"
        title = title or "Base frequency entropy"
        result = _plotBaseFrequenciesEntropy(
            significantOffsets,
            baseCountAtOffset,
            readCountAtOffset,
            outfile,
            title + subtitle,
            histogram,
            show,
            titleFontSize,
            axisFontSize,
        )

    if valuesFile:
        # The following will fail if sortOn is None (no result, above).
        with open(valuesFile, "w") as fp:
            dump(
                {
                    "sampleName": sampleName,
                    "text": [text for _, text in result],
                    "values": [value for value, _ in result],
                },
                fp,
            )


def plotCoverage(fig, row, col, readCountAtOffset, genomeLength):
    """
    Plot the read coverage along the genome.
    """
    meanCoverage = sum(readCountAtOffset) / genomeLength
    x = [i + 1 for i in range(genomeLength)]
    text = [str(i) for i in x]

    trace = go.Scatter(x=x, y=readCountAtOffset, showlegend=False, text=text)
    fig.append_trace(trace, row, col)

    # These are hacks. You shouldn't have to do things this way!
    fig["layout"]["annotations"][0]["text"] = (
        "Genome read coverage (mean %.3f)" % meanCoverage
    )
    fig["layout"]["yaxis1"].update({"title": "Read count"})
    fig["layout"]["xaxis"].update(
        {
            "range": (0, genomeLength + 1),
        }
    )
    fig["layout"]["yaxis"].update(
        {
            "range": (0, max(readCountAtOffset) + 1),
        }
    )


def plotSignificantOffsets(fig, row, col, significantOffsets, genomeLength):
    """
    Plot the genome offsets that are significant.
    """
    n = len(significantOffsets)
    trace = go.Scatter(
        x=[i + 1 for i in significantOffsets],
        y=[1.0] * n,
        mode="markers",
        showlegend=False,
    )
    fig.append_trace(trace, row, col)
    fig["layout"]["annotations"][1]["text"] = "%d significant genome location%s" % (
        n,
        s(n),
    )
    fig["layout"]["xaxis"].update(
        {
            "range": (0, genomeLength + 1),
        }
    )


def plotCoverageAndSignificantLocations(
    readCountAtOffset, genomeLength, significantOffsets, outfile, title=None, show=False
):
    """
    Plot read coverage and the significant locations.
    """
    fig = make_subplots(rows=2, cols=1, subplot_titles=("a", "b"), print_grid=False)

    plotCoverage(fig, 1, 1, readCountAtOffset, genomeLength)

    plotSignificantOffsets(fig, 2, 1, significantOffsets, genomeLength)

    if title is not None:
        fig["layout"].update(
            {
                "title": title,
            }
        )

    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def plotConsistentComponents(
    reference,
    components,
    outfile,
    infoFile,
    title="Consistent components",
    show=False,
    titleFontSize=12,
    axisFontSize=12,
    minReadsPerConsistentComponent=2,
):
    """
    Make a plot of all consistent connected components.

    @param reference: A C{Reference} instance.
    @param components: A C{list} of C{ComponentByOffsets} instances.
    @param outfile: The C{Path} to a file to write the plot to.
    @param infofile: The C{Path} to a file to write informative text output to.
    @param title: The C{str} title for the plot.
    @param titleFontSize: The C{int} title font size.
    @param axisFontSize: The C{int} axis font size.
    @param minReadsPerConsistentComponent: The C{int} minimum number of reads
        that must be in a consistent connected component in order for it to
        be shown.
    """

    def offsetsToLocationsStr(offsets):
        """
        Convert a list of zero-based offsets into a 1-based comma-separated string.

        @param offsets: An iterable of C{int} zero-based offsets.
        """
        return ", ".join(map(lambda i: str(i + 1), sorted(offsets)))

    data = []

    with open(infoFile, "w") as fp:
        print(
            "There are %d significant location%s: %s"
            % (
                len(reference.significantOffsets),
                s(len(reference.significantOffsets)),
                offsetsToLocationsStr(reference.significantOffsets),
            ),
            file=fp,
        )

        for count, component in enumerate(components, start=1):
            print(
                "Processing component %d, with %d consistent component%s"
                % (
                    count,
                    len(component.consistentComponents),
                    s(len(component.consistentComponents)),
                ),
                file=fp,
            )

            # Get the reference sequence for the component.
            reads = list(
                FastaReads(
                    reference.outputDir / ("component-%d-consensuses.fasta" % count)
                )
            )

            componentReferenceRead = reads[0]
            length = len(componentReferenceRead)
            minOffset = min(component.offsets)
            maxOffset = max(component.offsets)

            print(
                "  Offset range (inclusive): %d to %d" % (minOffset + 1, maxOffset),
                file=fp,
            )

            legendGroup = f"Component {count}"

            componentTotalReads = sum(
                len(cc.reads) for cc in component.consistentComponents
            )

            # Add a top line to represent the reference.
            data.append(
                go.Scatter(
                    x=(minOffset + 1, maxOffset + 1),
                    y=(1.05, 1.05),
                    hoverinfo="text",
                    name=legendGroup,
                    legendgroup=legendGroup,
                    text=(
                        "Component %d/%d: %d offset%s, %d reads, %d consistent "
                        "component%s"
                        % (
                            count,
                            len(components),
                            len(component.offsets),
                            s(len(component.offsets)),
                            componentTotalReads,
                            len(component.consistentComponents),
                            s(len(component.consistentComponents)),
                        )
                    ),
                )
            )

            # Add vertical lines at the start and end of this component.
            data.append(
                go.Scatter(
                    x=(minOffset + 1, minOffset + 1),
                    y=(-0.05, 1.05),
                    mode="lines",
                    hoverinfo="none",
                    line={
                        "color": "#ddd",
                    },
                    legendgroup=legendGroup,
                    showlegend=False,
                )
            )
            data.append(
                go.Scatter(
                    x=(maxOffset + 1, maxOffset + 1),
                    y=(-0.05, 1.05),
                    mode="lines",
                    hoverinfo="none",
                    line={
                        "color": "#ddd",
                    },
                    legendgroup=legendGroup,
                    showlegend=False,
                )
            )

            for ccCount, cc in enumerate(component.consistentComponents, start=1):
                if len(cc.reads) < minReadsPerConsistentComponent:
                    continue
                ccSummary = (
                    "Consistent component %d/%d. "
                    "Read count %d, offsets covered %d/%d"
                ) % (
                    ccCount,
                    len(component.consistentComponents),
                    len(cc.reads),
                    len(cc.nucleotides),
                    len(component.offsets),
                )

                # Get the consistent connected component consensus.
                consensus = reads[ccCount]
                assert ("consistent-component-%d" % ccCount) in consensus.id

                print("  Processing consistent component", ccCount, file=fp)
                print("  Component sequence:", consensus.sequence, file=fp)
                print(
                    "  %d offset%s: %s"
                    % (
                        len(cc.nucleotides),
                        s(len(cc.nucleotides)),
                        offsetsToLocationsStr(cc.nucleotides),
                    ),
                    file=fp,
                )

                match = compareDNAReads(componentReferenceRead, consensus)
                print(
                    matchToString(
                        match, componentReferenceRead, consensus, indent="    "
                    ),
                    file=fp,
                )

                identicalMatchCount = match["match"]["identicalMatchCount"]
                ambiguousMatchCount = match["match"]["ambiguousMatchCount"]

                # The match fraction will ignore gaps in the consensus
                # sequence as it is padded with '-' chars to align it to
                # the reference.
                fraction = (identicalMatchCount + ambiguousMatchCount) / (
                    length - len(match["read2"]["gapOffsets"])
                )

                x = []
                y = [fraction] * len(cc.nucleotides)
                text = []
                identical = []
                for index, offset in enumerate(sorted(component.offsets)):
                    if offset in cc.nucleotides:
                        consensusBase = consensus.sequence[index]
                        referenceBase = componentReferenceRead.sequence[index]

                        if consensusBase == referenceBase:
                            identical.append(len(x))

                        # x axis values are 1-based (locations, not offsets)
                        x.append(offset + 1)

                        text.append(
                            "%s<br>"
                            "Location: %d, component: %s, reference: %s<br>"
                            "Component nucleotides: %s"
                            % (
                                ccSummary,
                                offset + 1,
                                consensusBase,
                                referenceBase,
                                baseCountsToStr(cc.nucleotides[offset]),
                            )
                        )

                data.append(
                    go.Scatter(
                        x=x,
                        y=y,
                        hoverinfo="text",
                        selectedpoints=identical,
                        showlegend=False,
                        legendgroup=legendGroup,
                        text=text,
                        mode="markers",
                        selected={
                            "marker": {
                                "color": "green",
                                "size": 3,
                            }
                        },
                        unselected={
                            "marker": {
                                "color": "red",
                                "size": 3,
                            }
                        },
                    )
                )

    # Add the significant offsets.
    n = len(reference.significantOffsets)
    data.append(
        go.Scatter(
            x=[offset + 1 for offset in reference.significantOffsets],
            y=[-0.05] * n,
            text=[
                f"Significant site {offset + 1}"
                for offset in reference.significantOffsets
            ],
            hoverinfo="text",
            mode="markers",
            name=f"Significant sites ({n})",
        )
    )

    layout = go.Layout(
        title=title,
        titlefont={
            "size": titleFontSize,
        },
        xaxis={
            "range": (0, len(reference.read) + 1),
            "title": "Genome location",
            "titlefont": {
                "size": axisFontSize,
            },
            "showgrid": False,
        },
        yaxis={
            "range": (-0.1, 1.1),
            "title": "Nucleotide identity with reference sequence",
            "titlefont": {
                "size": axisFontSize,
            },
        },
        hovermode="closest",
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)

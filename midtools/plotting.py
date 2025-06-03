from __future__ import annotations

import sys
from pathlib import Path
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
from typing import Optional, TYPE_CHECKING

from dark.dna import compareDNAReads
from dark.fasta import FastaReads
from dark.sam import SAMFilter

if TYPE_CHECKING:
    from midtools.clusterAnalysis import Component
    from midtools.reference import Reference
from midtools.entropy import entropy2, MAX_ENTROPY
from midtools.match import matchToString
from midtools.utils import s, baseCountsToStr


try:
    # This relies on the private Cambridge/Charite HBV repo.
    from pyhbv.paperDieter import GENOTYPE_COLOR
    from pyhbv.genotype import getGenotype, genotypeKey
    from pyhbv.samples import sampleIdKey
except ImportError:
    haveHBV = False

    GENOTYPE_COLOR = {}

    def genotypeKey(genotype: str) -> tuple[int, str]:
        return 0, genotype

    def getGenotype(sample: str, trimSubgenotype: bool = True) -> Optional[str]:
        return None

    def sampleIdKey(id_: str) -> tuple[str, int, str]:
        return id_, 0, ""
else:
    haveHBV = True



def plotSAM(
    samFilter: SAMFilter,
    outfile: Path,
    title: str = "Reads",
    titleFontSize: int = 18,
    axisFontSize: int = 16,
    show: bool = False,
    jitter: float = 0.0,
) -> None:
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

    for alignment in samFilter.alignments():
        referenceStart = alignment.reference_start
        # Convert to float here because get_tag will return None if the tag
        # isn't present. This would cause the code to throw in any case, but
        # at least with a cast to float pyright doesn't complain.
        score = float(alignment.get_tag("AS")) + (
            0.0 if jitter == 0.0 else uniform(-jitter, jitter)
        )
        id_ = alignment.query_name
        assert isinstance(alignment.reference_length, int)
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
        "title": "Alignment score",
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
        hovermode="closest",
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def plotAllReferencesSAM(
    samFilter: SAMFilter,
    outfile: Path,
    sampleName: str,
    alignmentFile: Path,
    hbv: bool = False,
    jitter: float = 0.0,
    titleFontSize: int = 14,
    axisFontSize: int = 12,
    show: bool = False,
) -> None:
    """
    Plot the alignments found in a SAM file, even if the references are different
    lengths.

    @param sampleName: The C{str} name of the sample whose reads are being
        analysed.
    """
    if hbv and not haveHBV:
        raise ValueError(
            "You passed a true value for 'hbv' but the pyhbv module could not be "
            "imported."
        )

    referenceLengths = samFilter.referenceLengths()
    maxReferenceLength = max(referenceLengths.values())

    if hbv:
        # Ugly hack alert: If we are showing multiple HBV genomes, because
        # they have slightly different lengths I decided to scale X axis
        # positions slightly so that the (scaled) genomes have the same
        # length and when drawing horizontal lines for matching reads they
        # would be shown in the "same" place. With HBV genotypes, the maximum
        # difference between genome lengths is a touch over 2% (from 3182 to
        # 3248 nt), so this is a very minor visual adjustment. Maybe I
        # shouldn't have even done it.  And if we don't have HBV samples,
        # then the scale factor is 1.0 for all sequence ids (the keys of the
        # referenceLengths dict).
        referenceScaleFactor = {
            id_: maxReferenceLength / length
            for id_, length in referenceLengths.items()
        }
        # The referenceGenotype dict allows for accounting of HBV sample
        # genotypes and also the genotype labels in the legend.
        referenceGenotype = {
            id_: (getGenotype(id_) or "UNKNOWN") for id_ in referenceLengths
        }
    else:
        # Do no scaling of reference lengths.
        referenceScaleFactor = dict.fromkeys(referenceLengths, 1.0)
        # The referenceGenotype dict here just maps sample ids onto
        # themselves because we don't have actual genotypes. This variable
        # could have been better named. And it would have been better to
        # allow the caller to optionally pass us a sample-id to genotype
        # mapping.
        referenceGenotype = {id_: id_ for id_ in referenceLengths}

    genotypes = sorted(set(referenceGenotype.values()), key=genotypeKey)
    legendRank = {genotype: i for i, genotype in enumerate(genotypes)}
    nGenotypes = len(genotypes)
    genotypeReferences = defaultdict(set)

    data = []
    inLegend = set()
    genotypeCount: Counter[str] = Counter()

    readCount = 0
    for readCount, alignment in enumerate(samFilter.alignments(), start=1):
        referenceId = alignment.reference_name
        assert referenceId, f"No reference name was present in alignment {alignment!r}."
        scaleFactor = referenceScaleFactor[referenceId]
        start = alignment.reference_start
        assert isinstance(alignment.reference_length, int), (
            f"The reference length was not an integer in alignment {alignment!r}."
        )
        end = start + alignment.reference_length
        # Convert to float here because get_tag will return None if the tag
        # isn't present. This would cause the code to throw in any case, but
        # at least with a cast to float pyright doesn't complain.
        score = float(alignment.get_tag("AS")) + (
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
                for gt in sorted(genotypes)
                if genotypeCount[gt]
            )
            + "."
        )
        sampleGenotypeDesc = f"(genotype {getGenotype(sampleName)}) "
    else:
        genotypeReferencesDesc = sampleGenotypeDesc = ""

    title = "<br>".join(
        wrap(
            f"Best-matched genotypes for {readCount} reads for {sampleName} "
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

    if hbv:
        genotypeColor = GENOTYPE_COLOR
    else:
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
    significantOffsets: list[int],
    baseCountAtOffset: list[Counter[str]],
    readCountAtOffset: list[int],
    outfile: Path,
    title: str,
    histogram: bool,
    show: bool,
    titleFontSize: int,
    axisFontSize: int,
) -> list[tuple[float, str]]:
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
    significantOffsets: list[int],
    baseCountAtOffset: list[Counter[str]],
    outfile: Path,
    title: str,
    histogram: bool,
    show: bool,
    titleFontSize: int,
    axisFontSize: int,
) -> list[list[float], str]:
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
    genomeLength: int,
    significantOffsets: list[int],
    baseCountAtOffset: list[Counter[str]],
    readCountAtOffset: list[int],
    outfile: Path,
    title: str,
    show: bool,
    titleFontSize: int,
    axisFontSize: int,
    yRange: tuple[int, int],
) -> None:
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
    significantOffsets: list[int],
    baseCountAtOffset: list[Counter[str]],
    readCountAtOffset: list[int],
    outfile: Path,
    title: str,
    show: bool,
    titleFontSize: int,
    axisFontSize: int,
    yRange: tuple[float, float],
) -> None:
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
        title={
            "text": title,
            "font": {
                "size": titleFontSize,
            },
        },
        xaxis={
            "title": {
                "text": "Significant site index",
                "font": {
                    "size": axisFontSize,
                },
            },
        },
        yaxis={
            "title": {
                "text": "Nucleotide frequency",
                "font": {
                    "size": axisFontSize,
                },
            },
            "range": yRange,
        },
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=str(outfile), auto_open=show, show_link=False)


def plotBaseFrequencies(
    significantOffsets: list[int],
    baseCountAtOffset: list[Counter[str]],
    readCountAtOffset: list[int],
    outfile: Path,
    title: Optional[str] = None,
    titleFontSize: int = 12,
    axisFontSize: int = 12,
    yRange: tuple[float, float] = (0.0, 1.0),
    sampleName: Optional[str] = None,
    valuesFile: Optional[Path] = None,
    minReads: int = 5,
    homogeneousCutoff: float = 0.9,
    sortOn: Optional[str] = None,
    histogram: bool = False,
    show: bool = False,
) -> None:
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
        result = None
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
            outfile,
            title + subtitle,
            histogram,
            show,
            titleFontSize,
            axisFontSize,
        )

    if valuesFile:
        # Fail if there is no result (because sortOn is None, above).
        assert result is not None, (
            f"Logic error! sortOn in theory should be None, but is {sortOn!r}."
        )
        with open(valuesFile, "w") as fp:
            dump(
                {
                    "sampleName": sampleName,
                    "text": [text for _, text in result],
                    "values": [value for value, _ in result],
                },
                fp,
            )


def plotCoverage(
    fig, row: int, col: int, readCountAtOffset: list[int], genomeLength: int
) -> None:
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


def plotSignificantOffsets(
    fig, row: int, col: int, significantOffsets: list[int], genomeLength: int
) -> None:
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
    readCountAtOffset: list[int],
    genomeLength: int,
    significantOffsets: list[int],
    outfile: Path,
    title: Optional[str] = None,
    show: bool = False,
) -> None:
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
    reference: Reference,
    components: list[Component],
    outfile: Path,
    infoFile: Path,
    title: str = "Consistent components",
    show: bool = False,
    titleFontSize: int = 12,
    axisFontSize: int = 12,
    minReadsPerConsistentComponent: int = 2,
) -> None:
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

    def offsetsToLocationsStr(offsets: list[int]) -> str:
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
                    str(reference.outputDir / f"component-{count}-consensuses.fasta")
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
                        offsetsToLocationsStr(list(cc.nucleotides)),
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

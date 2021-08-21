from random import uniform
import plotly
from plotly import tools
from os.path import join
import plotly.graph_objs as go
from operator import itemgetter
from json import dump

from dark.dna import compareDNAReads
from dark.fasta import FastaReads

from midtools.entropy import entropy2, MAX_ENTROPY
from midtools.match import matchToString
from midtools.utils import s, baseCountsToStr


def plotSAM(samFilter, outfile, title='Reads', titleFontSize=18,
            axisFontSize=16, show=False, jitter=0.0):
    """
    Plot the alignments in a SAM file.
    """
    referenceLengths = samFilter.referenceLengths()

    if len(set(referenceLengths.values())) == 1:
        _, referenceLength = referenceLengths.popitem()
    else:
        raise ValueError(
            'SAM/BAM file reference sequences lengths (%s) are not '
            'all identical.' % ', '.join(map(str, sorted(referenceLengths))))

    data = []

    count = 0
    for count, alignment in enumerate(samFilter.alignments(), start=1):
        referenceStart = alignment.reference_start
        score = alignment.get_tag('AS') + (
            0.0 if jitter == 0.0 else uniform(-jitter, jitter))
        id_ = alignment.query_name
        data.append(go.Scatter(
            x=(referenceStart, referenceStart + alignment.query_length),
            y=(score, score),
            text=(id_, id_),
            hoverinfo='text',
            mode='lines',
            showlegend=False))

    xaxis = {
        'title': 'Genome location',
        'range': (0, referenceLength),
        'titlefont': {
            'size': axisFontSize,
        },
    }

    yaxis = {
        'title': 'score',
        # 'range': (0.0, 1.0),
        'titlefont': {
            'size': axisFontSize,
        },
    }

    title = '%s<br>%d reads%s' % (
        title, count,
        '' if jitter == 0.0 else (' (jitter %.2f)' % jitter))

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            'size': titleFontSize,
        },
        hovermode='closest',
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfile, auto_open=show, show_link=False)


def _plotSortedMaxBaseFrequencies(
        significantOffsets, baseCountAtOffset, readCountAtOffset, outfile,
        title, histogram, show, titleFontSize, axisFontSize):
    """
    Plot the sorted maximum base frequency for each of the significant
    offsets.
    """
    frequencyInfo = []

    for offset in significantOffsets:
        count = readCountAtOffset[offset]

        sortedFreqs = [x / count for x in
                       sorted(baseCountAtOffset[offset].values(),
                              reverse=True)]

        text = ('site %d<br>' % (offset + 1) +
                ', '.join('%s: %d' % (k, v)
                          for k, v in baseCountAtOffset[offset].items()))

        frequencyInfo.append((sortedFreqs[0], text))

    # We don't have to sort if we're making a histogram, but we're expected
    # to return a sorted values list, so we sort unconditionally.
    frequencyInfo.sort(key=itemgetter(0))
    values = [freq for freq, _ in frequencyInfo]

    if histogram:
        data = [
            go.Histogram(x=values, histnorm='probability'),
        ]

        xaxis = {
            'title': 'Significant site maximum nucleotide frequency',
            'range': (-0.05, 1.05),
            'titlefont': {
                'size': axisFontSize,
            },
        }

        yaxis = {
            'title': 'Probability mass',
            'range': (0.0, 1.0),
            'titlefont': {
                'size': axisFontSize,
            },
        }
    else:
        data = [
            go.Scatter(
                x=list(range(1, len(significantOffsets) + 1)),
                y=values,
                mode='markers',
                showlegend=False,
                text=[text for _, text in frequencyInfo]),
        ]

        xmargin = max(1, int(len(significantOffsets) * 0.01))
        xaxis = {
            'title': 'Rank',
            'range': (-xmargin, len(significantOffsets) + xmargin),
            'titlefont': {
                'size': axisFontSize,
            },
        }

        yaxis = {
            'range': (0.0, 1.05),
            'title': 'Significant site maximum nucleotide frequency',
            'titlefont': {
                'size': axisFontSize,
            },
        }

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            'size': titleFontSize,
        },
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfile, auto_open=show, show_link=False)
    return frequencyInfo


def _plotBaseFrequenciesEntropy(
        significantOffsets, baseCountAtOffset, readCountAtOffset, outfile,
        title, histogram, show, titleFontSize, axisFontSize):
    """
    Plot the sorted entropy of base frequencies for each of the significant
    offsets.
    """
    entropyInfo = []

    for offset in significantOffsets:
        text = ('site %d<br>' % (offset + 1) +
                ', '.join('%s: %d' % (k, v)
                          for k, v in baseCountAtOffset[offset].items()))

        entropyInfo.append(
            (entropy2(list(baseCountAtOffset[offset].elements())), text))

    assert all([ent <= MAX_ENTROPY for ent, _ in entropyInfo])

    # We don't have to sort if we're making a histogram, but we're expected
    # to return a sorted values list, so we sort unconditionally.
    entropyInfo.sort(key=itemgetter(0))
    values = [ent for ent, _ in entropyInfo]

    if histogram:
        data = [
            go.Histogram(x=values, histnorm='probability')
        ]

        xaxis = {
            'title': ('Significant site nucleotide frequency entropy '
                      '(bits)'),
            'range': (-0.05, MAX_ENTROPY),
            'titlefont': {
                'size': axisFontSize,
            },
        }

        yaxis = {
            'title': 'Probability mass',
            'range': (0.0, 1.0),
            'titlefont': {
                'size': axisFontSize,
            },
        }
    else:
        data = [
            go.Scatter(
                x=list(range(1, len(significantOffsets) + 1)),
                y=values,
                mode='markers',
                showlegend=False,
                text=[text for _, text in entropyInfo]),
        ]

        xmargin = max(1, int(len(significantOffsets) * 0.01))
        xaxis = {
            'range': (-xmargin, len(significantOffsets) + xmargin),
            'title': 'Rank',
            'titlefont': {
                'size': axisFontSize,
            },
        }

        yaxis = {
            'range': (-0.05, MAX_ENTROPY),
            'title': ('Significant site nucleotide frequency entropy '
                      '(bits)'),
            'titlefont': {
                'size': axisFontSize,
            },
        }

    layout = go.Layout(
        title=title,
        xaxis=xaxis,
        yaxis=yaxis,
        titlefont={
            'size': titleFontSize,
        },
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfile, auto_open=show, show_link=False)
    return entropyInfo


def _plotBaseFrequencies(significantOffsets, baseCountAtOffset,
                         readCountAtOffset, outfile, title, show,
                         titleFontSize, axisFontSize):
    """
    Plot the (sorted) base frequencies for each of the significant offsets.
    """
    x = list(range(len(significantOffsets)))
    text = []
    freqs = (
        [], [], [], [],
    )

    for offset in significantOffsets:
        count = readCountAtOffset[offset]

        sortedFreqs = [x / count for x in
                       sorted(baseCountAtOffset[offset].values(),
                              reverse=True)]
        while len(sortedFreqs) < 4:
            sortedFreqs.append(0.0)

        for i, frequency in enumerate(sortedFreqs):
            freqs[i].append(frequency)

        text.append(
            ('site %d<br>' % (offset + 1)) +
            ', '.join('%s: %d' % (k, v)
                      for k, v in baseCountAtOffset[offset].items()))

    data = [
        go.Bar(x=x, y=freqs[0], showlegend=False, text=text),
        go.Bar(x=x, y=freqs[1], showlegend=False),
        go.Bar(x=x, y=freqs[2], showlegend=False),
        go.Bar(x=x, y=freqs[3], showlegend=False),
    ]
    layout = go.Layout(
        barmode='stack',
        title=title,
        titlefont={
            'size': titleFontSize,
        },
        xaxis={
            'title': 'Significant site index',
            'titlefont': {
                'size': axisFontSize,
            },
        },
        yaxis={
            'title': 'Nucleotide frequency',
            'range': (0.45, 1.0),
            'titlefont': {
                'size': axisFontSize,
            },
        },
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfile, auto_open=show, show_link=False)


def plotBaseFrequencies(
        significantOffsets, baseCountAtOffset, readCountAtOffset, outfile,
        title=None, sampleName=None, valuesFile=None, minReads=5,
        homogeneousCutoff=0.9, sortOn=None, histogram=False, show=False,
        titleFontSize=12, axisFontSize=12):
    """
    Plot sorted base frequencies at signifcant sites.
    """

    subtitle = (
        '<br>%d significant sites. Min %d read%s per site. '
        '%.2f homogeneity cutoff.' %
        (len(significantOffsets), minReads, s(minReads), homogeneousCutoff))

    if sortOn is None:
        title = title or 'Base frequencies (sorted)'
        _plotBaseFrequencies(significantOffsets, baseCountAtOffset,
                             readCountAtOffset, outfile,
                             title + subtitle, show,
                             titleFontSize, axisFontSize)
    elif sortOn == 'max':
        title = title or 'Maximum base frequency'
        result = _plotSortedMaxBaseFrequencies(
            significantOffsets, baseCountAtOffset,
            readCountAtOffset, outfile, title + subtitle,
            histogram, show, titleFontSize, axisFontSize)
    else:
        assert sortOn == 'entropy', (
            'Unknown --sortOn value: %r' % sortOn)
        title = title or 'Base frequency entropy'
        result = _plotBaseFrequenciesEntropy(
            significantOffsets, baseCountAtOffset,
            readCountAtOffset, outfile, title + subtitle,
            histogram, show, titleFontSize, axisFontSize)

    if valuesFile:
        # The following will fail if sortOn is None (no result, above).
        with open(valuesFile, 'w') as fp:
            dump(
                {
                    'sampleName': sampleName,
                    'text': [text for _, text in result],
                    'values': [value for value, _ in result],
                },
                fp
            )


def plotCoverage(fig, row, col, readCountAtOffset, genomeLength):
    """
    Plot the read coverage along the genome.
    """
    meanCoverage = sum(readCountAtOffset) / genomeLength
    x = [i + 1 for i in range(genomeLength)]
    text = [str(i) for i in x]

    trace = go.Scatter(
        x=x,
        y=readCountAtOffset,
        showlegend=False,
        text=text)
    fig.append_trace(trace, row, col)

    # These are hacks. You shouldn't have to do things this way!
    fig['layout']['annotations'][0]['text'] = (
        'Genome read coverage (mean %.3f)' % meanCoverage)
    fig['layout']['yaxis1'].update({
        'title': 'Read count'
    })
    fig['layout']['xaxis'].update({
        'range': (0, genomeLength + 1),
    })
    fig['layout']['yaxis'].update({
        'range': (0, max(readCountAtOffset) + 1),
    })


def plotSignificantOffsets(fig, row, col, significantOffsets, genomeLength):
    """
    Plot the genome offsets that are significant.
    """
    n = len(significantOffsets)
    trace = go.Scatter(
        x=[i + 1 for i in significantOffsets],
        y=[1.0] * n,
        mode='markers',
        showlegend=False)
    fig.append_trace(trace, row, col)
    fig['layout']['annotations'][1]['text'] = (
        '%d significant genome location%s' % (n, s(n)))
    fig['layout']['xaxis'].update({
        'range': (0, genomeLength + 1),
    })


def plotCoverageAndSignificantLocations(
        readCountAtOffset, genomeLength, significantOffsets, outfile,
        title=None, show=False):
    """
    Plot read coverage and the significant locations.
    """
    fig = tools.plotly.subplots.make_subplots(
        rows=2, cols=1, subplot_titles=('a', 'b'), print_grid=False)

    plotCoverage(fig, 1, 1, readCountAtOffset, genomeLength)

    plotSignificantOffsets(fig, 2, 1, significantOffsets, genomeLength)

    if title is not None:
        fig['layout'].update({
            'title': title,
        })

    plotly.offline.plot(fig, filename=outfile, auto_open=show,
                        show_link=False)


def plotConsistentComponents(
        referenceId, genomeLength, components, significantOffsets, outfile,
        infoFile, outputDir, title='xxx', show=False, titleFontSize=12,
        axisFontSize=12):
    """
    Plot consistent connected components.
    """
    def offsetsToLocationsStr(offsets):
        return ', '.join(map(lambda i: str(i + 1), sorted(offsets)))

    data = []

    with open(infoFile, 'w') as fp:

        print('There are %d significant location%s: %s' %
              (len(significantOffsets), s(len(significantOffsets)),
               offsetsToLocationsStr(significantOffsets)), file=fp)

        for count, component in enumerate(components, start=1):

            print('Processing component %d, with %d consistent component%s' %
                  (count, len(component.consistentComponents),
                   s(len(component.consistentComponents))), file=fp)

            # Get the reference sequence for the component.
            reads = list(FastaReads(
                join(outputDir, 'component-%d-consensuses.fasta' % count)))

            reference = reads[0]
            length = len(reference)
            minOffset = min(component.offsets)
            maxOffset = max(component.offsets)

            print('  Offset range: %d to %d' % (minOffset + 1, maxOffset + 1),
                  file=fp)

            # Add a top line to represent the reference.
            data.append(go.Scatter(
                x=(minOffset + 1, maxOffset + 1),
                y=(1.05, 1.05),
                hoverinfo='text',
                name=('Reference component %s' % count),
                text=('Reference component %s, %d offsets' %
                      (count, len(component.offsets)))))

            # Add vertical lines at the start and end of this component.
            data.append(go.Scatter(
                x=(minOffset + 1, minOffset + 1),
                y=(-0.05, 1.05),
                mode='lines',
                hoverinfo='none',
                line={
                    'color': '#eee',
                },
                showlegend=False,
            ))
            data.append(go.Scatter(
                x=(maxOffset + 1, maxOffset + 1),
                y=(-0.05, 1.05),
                mode='lines',
                hoverinfo='none',
                line={
                    'color': '#eee',
                },
                showlegend=False,
            ))

            for ccCount, cc in enumerate(component.consistentComponents,
                                         start=1):

                ccSummary = ('Component read count %d, offsets covered %d/%d' %
                             (len(cc.reads), len(cc.nucleotides),
                              len(component.offsets)))

                # Get the consistent connected component consensus.
                consensus = reads[ccCount]
                assert ('consistent-component-%d' % ccCount) in consensus.id

                print('  Processing consistent component', ccCount, file=fp)
                print('  Component sequence:', consensus.sequence, file=fp)
                print('  %d offset%s: %s' %
                      (len(cc.nucleotides), s(len(cc.nucleotides)),
                       offsetsToLocationsStr(cc.nucleotides)), file=fp)

                match = compareDNAReads(reference, consensus)
                print(
                    matchToString(match, reference, consensus, indent='    '),
                    file=fp)

                identicalMatchCount = match['match']['identicalMatchCount']
                ambiguousMatchCount = match['match']['ambiguousMatchCount']

                # The match fraction will ignore gaps in the consensus
                # sequence as it is padded with '-' chars to align it to
                # the reference.
                fraction = (identicalMatchCount + ambiguousMatchCount) / (
                    length - len(match['read2']['gapOffsets']))

                x = []
                y = [fraction] * len(cc.nucleotides)
                text = []
                identical = []
                for index, offset in enumerate(sorted(component.offsets)):
                    if offset in cc.nucleotides:

                        consensusBase = consensus.sequence[index]
                        referenceBase = reference.sequence[index]

                        if consensusBase == referenceBase:
                            identical.append(len(x))

                        # x axis values are 1-based (locations, not offsets)
                        x.append(offset + 1)

                        text.append(
                            'Location: %d, component: %s, reference: %s'
                            '<br>Component nucleotides: %s<br>%s' %
                            (offset + 1,
                             consensusBase,
                             referenceBase,
                             baseCountsToStr(cc.nucleotides[offset]),
                             ccSummary))

                data.append(
                    go.Scatter(
                        x=x,
                        y=y,
                        hoverinfo='text',
                        selectedpoints=identical,
                        showlegend=False,
                        text=text,
                        mode='markers',
                        selected={
                            'marker': {
                                'color': 'blue',
                            }
                        },
                        unselected={
                            'marker': {
                                'color': 'red',
                            }
                        }
                    )
                )

    # Add the significant offsets.
    n = len(significantOffsets)
    data.append(go.Scatter(
        x=[i + 1 for i in significantOffsets],
        y=[-0.05] * n,
        text=['Location %d' % (offset + 1) for offset in significantOffsets],
        hoverinfo='text',
        mode='markers',
        name='Significant locations'))

    layout = go.Layout(
        title=title,
        titlefont={
            'size': titleFontSize,
        },
        xaxis={
            'range': (0, genomeLength + 1),
            'title': 'Genome location',
            'titlefont': {
                'size': axisFontSize,
            },
        },
        yaxis={
            'range': (-0.1, 1.1),
            'title': 'Nucleotide identity with reference sequence',
            'titlefont': {
                'size': axisFontSize,
            },
        },
        hovermode='closest',
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfile, auto_open=show, show_link=False)

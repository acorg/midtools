import plotly
from plotly import tools
import plotly.graph_objs as go
from operator import itemgetter
from json import dump

from mid.entropy import entropy2, MAX_ENTROPY
from mid.utils import s


def _plotSortedMaxBaseFrequencies(
        significantOffsets, baseCountAtOffset, readCountAtOffset, outFile,
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
    plotly.offline.plot(fig, filename=outFile, auto_open=show, show_link=False)
    return frequencyInfo


def _plotBaseFrequenciesEntropy(
        significantOffsets, baseCountAtOffset, readCountAtOffset, outFile,
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
    plotly.offline.plot(fig, filename=outFile, auto_open=show, show_link=False)
    return entropyInfo


def _plotBaseFrequencies(significantOffsets, baseCountAtOffset,
                         readCountAtOffset, outFile, title, show,
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
    plotly.offline.plot(fig, filename=outFile, auto_open=show, show_link=False)


def plotBaseFrequencies(
        significantOffsets, baseCountAtOffset, readCountAtOffset, outFile,
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
                             readCountAtOffset, outFile,
                             title + subtitle, show,
                             titleFontSize, axisFontSize)
    elif sortOn == 'max':
        title = title or 'Maximum base frequency'
        result = _plotSortedMaxBaseFrequencies(
            significantOffsets, baseCountAtOffset,
            readCountAtOffset, outFile, title + subtitle,
            histogram, show, titleFontSize, axisFontSize)
    else:
        assert sortOn == 'entropy', (
            'Unknown --sortOn value: %r' % sortOn)
        title = title or 'Base frequency entropy'
        result = _plotBaseFrequenciesEntropy(
            significantOffsets, baseCountAtOffset,
            readCountAtOffset, outFile, title + subtitle,
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
        x=x, y=readCountAtOffset, showlegend=False, text=text)
    fig.append_trace(trace, row, col)

    # These are hacks. You shouldn't have to do things this way!
    fig['layout']['annotations'][0]['text'] = (
        'Genome read coverage (mean %.3f)' % meanCoverage)
    fig['layout']['yaxis1'].update({
        'title': 'Read count'
    })


def plotSignificantOffsets(fig, row, col, significantOffsets):
    """
    Plot the genome offsets that are significant.
    """
    n = len(significantOffsets)
    trace = go.Scatter(
        x=[i + 1 for i in significantOffsets], y=[1.0] * n,
        mode='markers', showlegend=False)
    fig.append_trace(trace, row, col)
    fig['layout']['annotations'][1]['text'] = (
        '%d significant genome location%s' % (n, s(n)))


def plotCoverageAndSignificantLocations(
        readCountAtOffset, genomeLength, significantOffsets, outFile,
        title=None, show=False):
    """
    Plot read coverage and the significant locations.
    """
    fig = tools.make_subplots(rows=2, cols=1, subplot_titles=('a', 'b'),
                              print_grid=False)

    plotCoverage(fig, 1, 1, readCountAtOffset, genomeLength)

    plotSignificantOffsets(fig, 2, 1, significantOffsets)

    if title is not None:
        fig['layout'].update({
            'title': title,
        })

    plotly.offline.plot(fig, filename=outFile, auto_open=show,
                        show_link=False)

#!/usr/bin/env python

from __future__ import division, print_function

import plotly
import plotly.graph_objs as go
from plotly import tools

from data.data import (
    addCommandLineOptions, parseCommandLineOptions, findSignificantOffsets)


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
        '%d significant genome location%s' % (n, '' if n == 1 else 's'))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Analyze a set of aligned reads.')

    parser.add_argument(
        '--title',
        help='The overall title plot')

    addCommandLineOptions(parser, 'coverage-and-significant-locations.html')
    args = parser.parse_args()

    (genome, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset) = parseCommandLineOptions(args)

    significantOffsets = list(findSignificantOffsets(
        baseCountAtOffset, readCountAtOffset, args.minReads,
        args.homogeneousCutoff))

    print('Read %d aligned reads. Found %d significant locations.' %
          (len(alignedReads), len(significantOffsets)))

    fig = tools.make_subplots(rows=2, cols=1,
                              subplot_titles=('a', 'b'))

    plotCoverage(fig, 1, 1, readCountAtOffset, len(genome))

    plotSignificantOffsets(fig, 2, 1, significantOffsets)

    if args.title is not None:
        fig['layout'].update({
            'title': args.title,
        })

    plotly.offline.plot(fig, filename=args.outFile, auto_open=args.show,
                        show_link=False)

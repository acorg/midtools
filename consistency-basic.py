#!/usr/bin/env python

from __future__ import division, print_function

import plotly
import plotly.graph_objs as go

from random import uniform
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score

from data.data import addCommandLineOptions, parseCommandLineOptions


def plotConsistency(significantOffsets, baseCountAtOffset,
                    readsAtOffset, outFile, title, jitter, show):
    """
    """
    x = []
    ariScores = []
    amiScores = []
    text = []

    for index in range(len(significantOffsets) - 1):
        offset1, offset2 = significantOffsets[index:index + 2]
        readsCoveringBoth = (readsAtOffset[offset1] &
                             readsAtOffset[offset2])
        # print('analyzeConsistency', offset1, offset2)
        # print('%d reads cover 1, %d cover 2, %d cover both' % (
        #     len(readsAtOffset[offset1]),
        #     len(readsAtOffset[offset2]),
        #     len(readsCoveringBoth)))

        if readsCoveringBoth:
            bases1 = []
            bases2 = []
            for read in readsCoveringBoth:
                bases1.append(read.base(offset1))
                bases2.append(read.base(offset2))

            x.append(offset1 + 1)
            ariScores.append(
                adjusted_rand_score(bases1, bases2) +
                (uniform(-0.01, 0.01) if jitter else 0.0)
            )
            text.append('location=%d' % x[-1])
            amiScores.append(
                adjusted_mutual_info_score(bases1, bases2) +
                (uniform(-0.01, 0.01) if jitter else 0.0)
            )

    data = [
        go.Scatter(x=x, y=ariScores, mode='markers', text=text, name='ARI'),
        go.Scatter(x=x, y=amiScores, mode='markers', text=text, name='AMI')
    ]

    if title is None:
        layout = go.Layout()
    else:
        layout = go.Layout({
            'title': title,
        })

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outFile, auto_open=show, show_link=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Analyze a set of aligned reads.')

    parser.add_argument(
        '--title',
        help='The title plot')

    parser.add_argument(
        '--jitter', action='store_true', default=False,
        help='Add jitter to scores to avoid over-plotting identical values.')

    addCommandLineOptions(parser, 'consistency-basic.html')
    args = parser.parse_args()

    (genomeLength, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset,
     significantOffsets) = parseCommandLineOptions(args, True)

    print('Read %d aligned reads. Found %d significant locations.' %
          (len(alignedReads), len(significantOffsets)))

    title = args.title or (
        'ARI and AMI at %d significant locations (%s).' %
        (len(significantOffsets),
         'jitter added' if args.jitter else 'no jitter'))

    plotConsistency(significantOffsets, baseCountAtOffset, readsAtOffset,
                    args.outFile, title, args.jitter, args.show)

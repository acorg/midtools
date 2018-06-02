#!/usr/bin/env python

from __future__ import division, print_function

import plotly
import plotly.graph_objs as go

from sklearn.metrics import adjusted_rand_score  # , adjusted_mutual_info_score

from data.data import addCommandLineOptions, parseCommandLineOptions
from data.nid import normalized_information_distance


def zeroScore(_, __):
    return 0.0


def plotConsistency(significantOffsets, baseCountAtOffset,
                    readsAtOffset, minCommonReads, outFile, title, show):
    """
    """
    scores = []
    text = []

    for offset1 in significantOffsets:
        rowScores = []
        rowText = []
        for offset2 in significantOffsets:
            if offset2 == offset1:
                rowScores.append(1.0)
                rowText.append(
                    '%d identity (%d read%s)' %
                    (offset1 + 1,
                        len(readsAtOffset[offset1]),
                     '' if len(readsAtOffset[offset1]) == 1 else 's'))
            else:
                readsCoveringBoth = (readsAtOffset[offset1] &
                                     readsAtOffset[offset2])

                # print('analyzeConsistency', offset1, offset2)
                # print('%d reads cover 1, %d cover 2, %d cover both' % (
                #     len(readsAtOffset[offset1]),
                #     len(readsAtOffset[offset2]),
                #     len(readsCoveringBoth)))

                if readsCoveringBoth:
                    commonCount = len(readsCoveringBoth)
                    if commonCount < minCommonReads:
                        rowScores.append(0.0)
                        rowText.append(
                            '%d (%d read%s) vs %d (%d read%s)<br>Too few (%d) reads in common' %
                            (offset2 + 1, len(readsAtOffset[offset2]), '' if len(readsAtOffset[offset2]) == 1 else 's',
                             offset1 + 1, len(readsAtOffset[offset1]), '' if len(readsAtOffset[offset1]) == 1 else 's',
                            commonCount))
                    else:
                        bases1 = []
                        bases2 = []
                        for read in readsCoveringBoth:
                            bases1.append(read.base(offset1))
                            bases2.append(read.base(offset2))

                        if offset1 == 1895 and offset2 == 1816:
                            print('offset1 == 1896 and offset2 == 1817')
                            print(bases1)
                            print(bases2)
                        # f = (adjusted_mutual_info_score if offset2 > offset1
                        #      else adjusted_rand_score)
                        f = (normalized_information_distance if offset2 > offset1
                             else adjusted_rand_score)

                        rowScores.append(f(bases1, bases2))
                        rowText.append(
                            '%d (%d read%s) vs %d (%d read%s)<br>%d read%s in common' %
                            (offset2 + 1, len(readsAtOffset[offset2]), '' if len(readsAtOffset[offset2]) == 1 else 's',
                             offset1 + 1, len(readsAtOffset[offset1]), '' if len(readsAtOffset[offset1]) == 1 else 's',
                             commonCount,
                             '' if commonCount == 1 else 's'))
                        # print('offset %d vs %d = %.4f (%d)' % (offset1, offset2, rowScores[-1], commonCount))
                else:
                    rowScores.append(-0.25)
                    rowText.append(
                        '%d vs %d, no reads in common' % (offset1 + 1, offset2 + 1))

        scores.append(rowScores)
        text.append(rowText)

    data = [
        # go.Heatmap(x=x, y=x, z=scores, name='ARI'),
        go.Heatmap(z=scores, name='ARI', text=text)  # , colorscale='Viridis'),
    ]

    layoutDict = dict(
        xaxis={
            'title': 'Significant location index',
        },
        yaxis={
            'title': 'Score',
        },
    )

    if title:
        layoutDict['title'] = title

    layout = go.Layout(layoutDict)

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outFile, auto_open=show, show_link=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Analyze a set of aligned reads.')

    addCommandLineOptions(parser, 'consistency-heatmap.html')

    parser.add_argument(
        '--title',
        help='The title plot')

    parser.add_argument(
        '--minCommonReads', type=int, default=5,
        help=('The minimum number of reads that a pair of locations must '
              'share for them to be considered in the consistency '
              'calculation.'))

    args = parser.parse_args()

    (genomeLength, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset,
     significantOffsets) = parseCommandLineOptions(args, True)

    print('Read %d aligned reads of length %d. '
          'Found %d significant offsets.' %
          (len(alignedReads), genomeLength, len(significantOffsets)))

    title = args.title or (
        'Consistency heatmap at %d significant offsets.<br>ARI above '
        'diagonal, NID below' % len(significantOffsets))

    plotConsistency(significantOffsets, baseCountAtOffset, readsAtOffset,
                    args.minCommonReads, args.outFile, title, args.show)

#!/usr/bin/env python

from __future__ import division, print_function

from os.path import basename
import plotly
import plotly.graph_objs as go

from dark.reads import readClassNameToClass

from data.data import (
    AlignedRead, addCommandLineOptions, parseCommandLineOptions,
    findSignificantOffsets, gatherData)


def baseCountsToStr(counts):
    """
    @param counts: A C{counter} instance.
    """
    return ' '.join([
        ('%s:%d' % (base, counts[base])) for base in sorted(counts)])


def plotBaseFrequencies(significantOffsets, baseCountAtOffset,
                        readCountAtOffset, outFile, title):
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
            ('location %d<br>' % (offset + 1)) +
            ', '.join('%s: %d' % (k, v)
                      for k, v in baseCountAtOffset[offset].items()))

        # print('loc %d, counts %s, freqs %s' % (offset,
        # baseCountAtOffset[offset], sortedFreqs))

    data = [
        go.Bar(x=x, y=freqs[0], showlegend=False, text=text),
        go.Bar(x=x, y=freqs[1], showlegend=False),
        go.Bar(x=x, y=freqs[2], showlegend=False),
        go.Bar(x=x, y=freqs[3], showlegend=False),
    ]
    layout = go.Layout(
        barmode='stack',
        title=title,
        xaxis={
            'title': 'Significant location index',
        },
        yaxis={
            'title': 'Nucleotide frequency',
        },
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outFile, show_link=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Analyze a set of aligned reads.')

    addCommandLineOptions(parser, 'significant-base-frequencies.html')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print verbose textual output showing read connections.')

    parser.add_argument(
        '--alignmentFasta',
        help=('The filename of an optional FASTA alignment to show the base '
              'frequencies for. The significant locations will always be '
              'calculated from the aligned reads, but you can use this '
              'argument to show the frequencies from another alignment at '
              'those locations.'))

    args = parser.parse_args()

    (genome, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset) = parseCommandLineOptions(args)

    genomeLength = len(genome)

    significantOffsets = list(findSignificantOffsets(
        baseCountAtOffset, readCountAtOffset, args.minReads,
        args.homogeneousCutoff))

    print('Read %d aligned reads of length %d. '
          'Found %d significant locations.' %
          (len(alignedReads), genomeLength, len(significantOffsets)))

    if args.alignmentFasta:
        # Now we have the significant offsets from the reads aligned to
        # the consensus, read in another alignment so we can show its base
        # frequencies at those offsets.
        readClass = readClassNameToClass[args.readClass]
        if args.fasta:
            from dark.fasta import FastaReads
            reads = FastaReads(args.alignmentFasta, readClass=readClass)
        elif args.fastq:
            from dark.fastq import FastqReads
            reads = FastqReads(args.alignmentFasta, readClass=readClass)
        else:
            from dark.fasta_ss import SSFastaReads
            reads = SSFastaReads(args.alignmentFasta, readClass=readClass)

        alignedReads = [AlignedRead(read) for read in reads
                        if len(read) == genomeLength]
        print('Read %d sequences of length %d from %s.' % (
            len(alignedReads), genomeLength, args.alignmentFasta))

        altReadCountAtOffset, altBaseCountAtOffset, altReadsAtOffset = (
            gatherData(genomeLength, alignedReads))

        if args.verbose:
            # Print a comparitive summary of bases at the significant
            # offsets.
            print('Summary of significant location base frequencies:')
            for offset in significantOffsets:
                print('  ', offset + 1, 'ref:',
                      baseCountsToStr(altBaseCountAtOffset[offset]),
                      'reads:',
                      baseCountsToStr(baseCountAtOffset[offset]))

        readCountAtOffset, baseCountAtOffset, readsAtOffset = (
            altReadCountAtOffset, altBaseCountAtOffset, altReadsAtOffset)

        title = ('Base frequencies of %d sequences from %s' %
                 (len(alignedReads), basename(args.alignmentFasta)))
    else:
        title = 'Base frequencies of %d aligned ancient read%s' % (
            len(alignedReads), '' if len(alignedReads) == 1 else 's')

    title += ('<br>at %d significant locations. Min %d read%s per location.' %
              (len(significantOffsets), args.minReads,
               '' if args.minReads == 1 else 's'))

    if args.trim:
        title += ' Reads trimmed by %d base%s.' % (
            args.trim, '' if args.trim == 1 else 's')

    if args.show:
        plotBaseFrequencies(significantOffsets, baseCountAtOffset,
                            readCountAtOffset, args.outFile, title)

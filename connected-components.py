#!/usr/bin/env python

from __future__ import division, print_function

import sys
import plotly
import plotly.graph_objs as go
import re
import colorlover
from collections import Counter
from pathlib import Path  # This is Python 3 only.
from os.path import join
from math import log10
from os import unlink
from itertools import chain

from collections import defaultdict
# from time import time, ctime
# import sys

from data.data import (
    addCommandLineOptions, parseCommandLineOptions, findSignificantOffsets)
from data.graph import Node, connectedComponents, componentOffsets
from dark.reads import Read

# The following does not include a code for all possible sets of
# nucleotides. Why?
AMBIGUOUS = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'M': {'A', 'C'},
    'R': {'A', 'G'},
    'W': {'A', 'T'},
    'S': {'G', 'C'},
    'K': {'G', 'T'},
    'Y': {'C', 'T'},
    'V': {'A', 'C', 'G'},
    'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'},
    'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'},
}

# Make a reverse version of AMBIGUOUS.
BASES_TO_AMBIGUOUS = {}
for symbol, bases in AMBIGUOUS.items():
    BASES_TO_AMBIGUOUS[''.join(sorted(bases))] = symbol


def baseCountsToStr(counts):
    """
    @param counts: A C{counter} instance.
    """
    return ' '.join([
        ('%s:%d' % (base, counts[base])) for base in sorted(counts)])


def plotComponents(components, categoryRegexs, categoryNames, genomeLength,
                   significantOffsets, outFile, title, additionalOffsets,
                   additionalOffsetsTitle):
    """
    Show all connected components in a plot, as well as the significant
    locations and any additional locations.
    """
    offsetsForLevel = []
    compiledCategoryRegexs = list(map(re.compile, categoryRegexs))
    nCategories = len(categoryRegexs)
    nColorCategories = max(
        3,
        1 +  # Significant offsets.
        (1 if additionalOffsets else 0) +  # Additional offsets.
        len(categoryRegexs))  # Read categories.
    colors = colorlover.scales[str(nColorCategories)]['qual']['Paired']
    shownLegendColorIndices = set()

    # Show the significant offsets.
    data = [
        go.Scatter(
            x=list(map(lambda offset: offset + 1, significantOffsets)),
            y=[-0.1] * len(significantOffsets),
            mode='markers',
            hoverinfo='x',
            marker={'color': colors[-1]},
            name='Significant location')
    ]

    # Show the additional offsets, if any.
    if additionalOffsets:
        data.append(go.Scatter(
            x=list(map(lambda offset: offset + 1, additionalOffsets)),
            y=[-0.2] * len(fields),
            mode='markers',
            hoverinfo='x',
            marker={'color': colors[-2]},
            name=additionalOffsetsTitle))

    for componentIndex, (minOffset, maxOffset, component) in enumerate(
            components):
        componentSize = len(component)
        # print('component %d len %d' % (componentIndex + 1, componentSize))
        # print('  minOffset %d, maxOffset %d, range %d' % (
        # minOffset, maxOffset, maxOffset - minOffset))
        offsets = set(range(minOffset, maxOffset))
        level = 0
        for level, levelOffsets in enumerate(offsetsForLevel):
            if not (levelOffsets & offsets):
                # This component's offsets can all be shown on this level.
                levelOffsets.update(offsets)
                break
        else:
            # Didn't break, so it's a new level.
            offsetsForLevel.append(offsets)
            level = len(offsetsForLevel) - 1

        # One extra slot in the following to count reads whose ids don't
        # match any category regex.
        categoryCounts = [0] * (nCategories + 1)

        for node in component:
            id_ = node.read.read.id
            for i, regex in enumerate(compiledCategoryRegexs):
                if regex.match(id_):
                    categoryCounts[i] += 1
                    break
            else:
                # No matching category for this read id!
                categoryCounts[-1] += 1

        # Sanity check that every read is in a category.
        assert sum(categoryCounts) == componentSize

        categorySummary = []
        for categoryIndex, categoryCount in enumerate(categoryCounts):
            if categoryCount:
                categorySummary.append(
                    'Category: %s. %d read%s' % (
                        categoryNames[categoryIndex], categoryCount,
                        '' if categoryCount == 1 else 's'))
        categorySummary = '<br>'.join(categorySummary)

        scale = (maxOffset - minOffset) / componentSize
        # print('  scale is %f' % scale)
        # print('  categoryCounts', categoryCounts)
        cumulative = 0

        for colorIndex, categoryCount in enumerate(categoryCounts):
            if categoryCount:
                increment = int(categoryCount * scale + 0.5)
                # print('    increment', increment)
                x = [minOffset + cumulative + 1,
                     minOffset + cumulative + increment + 1]
                y = [level, level]
                # print('   ', x)
                cumulative += increment
                data.append(go.Scatter(
                    x=x, y=y, mode='lines',
                    showlegend=(colorIndex not in shownLegendColorIndices),
                    legendgroup=str(colorIndex),
                    text=('Component %d, %d read%s<br>Len %d (locations %d-%d)'
                          '<br>%s' % (
                              componentIndex + 1, componentSize,
                              '' if componentSize == 1 else 's',
                              maxOffset - minOffset,
                              minOffset + 1, maxOffset, categorySummary)),
                    hoveron='points+fills',
                    name=categoryNames[colorIndex],
                    hoverinfo='text',
                    line={
                        'color': colors[colorIndex],
                        'width': int(log10(componentSize) + 1.0) * 10,
                    }))
                shownLegendColorIndices.add(colorIndex)

    layout = go.Layout(
        xaxis={
            'visible': True,
            'title': 'Genome position',
            'range': (0, genomeLength),
        },
        yaxis={
            'visible': False,
            'zeroline': False,
            'showline': False,
            'ticks': '',
            'showticklabels': False,
            'autorange': True,
            'showgrid': False,
        },
        hovermode='closest',
    )
    if title:
        layout['title'] = title

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outFile, show_link=False)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Find which reads agree and disagree with one another.')

    addCommandLineOptions(parser, 'connected-components.html')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print verbose textual output showing read connections.')

    parser.add_argument(
        '--componentJSON',
        help='The filename to write JSON containing all components to.')

    parser.add_argument(
        '--fastaDir',
        help=('The directory to write FASTA of all components to. The files '
              'will be named component-1.fasta, component-2.fasta, etc. '
              'All pre-existing files with such names *will be removed* to '
              'avoid confusion due to left-over files from earlier runs. '
              'If not specified, component FASTA files will not be written.'))

    parser.add_argument(
        '--additionalLocations',
        help=('Specify some additional locations to highlight. The argument '
              'must be a string name of the locations followed by space '
              'separated genome locations (1-based). Underscores in the '
              'name will be replaced by spaces.'))

    parser.add_argument(
        '--readCategoryRegex', nargs='*',
        help=('Specify regular expressions that match read ids for the '
              'purpose of plotting components and showing their catgeory '
              'composition. The category will normally be something that '
              'indicates to which genome a read belongs.'))

    parser.add_argument(
        '--readCategoryRegexNames', nargs='*',
        help=('Specify names for the read categories defined in the regular '
              'expressions in --readCategoryRegex. If not specified, the '
              'regular expressions (i.e., their strings) will be used as the '
              'category names. If specified, must be the same length as '
              '--readCategoryRegex.'))

    parser.add_argument(
        '--agreementCutoff', type=float, default=0.9,
        help=('Only reads that agree at least this much will be considered '
              'connected.'))

    parser.add_argument(
        '--runSpades', action='store_true', default=False,
        help=('If specified, run the SPAdes assembler on the reads in each '
              'component. The output directories will be named '
              'component-1-spades etc., and will be created in the directory '
              'specified by --fastaDir (which must therefore also be '
              'specified).'))

    parser.add_argument(
        '--addConsensusToAlignment', action='store_true', default=False,
        help=('If specified, put the consensus sequence into the alignment '
              'for each connected component (only applies if FASTA is '
              'being written for each component, dut to --fastaDir being '
              'specified).'))

    args = parser.parse_args()

    if args.readCategoryRegexNames:
        if len(args.readCategoryRegexNames) != len(args.readCategoryRegex):
            print('--readCategoryRegex and --readCategoryRegexNames must be '
                  'of equal length.', file=sys.stderr)
            sys.exit(1)

        readCategoryRegexNames = args.readCategoryRegexNames
    else:
        readCategoryRegexNames = list(args.readCategoryRegex)

    if args.runSpades and not args.fastaDir:
        print('--fastaDir must be specified if --runSpades is used.',
              file=sys.stderr)
        sys.exit(1)

    (genome, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset) = parseCommandLineOptions(args)

    genomeLength = len(genome)

    significantOffsets = list(findSignificantOffsets(
        baseCountAtOffset, readCountAtOffset, args.minReads,
        args.homogeneousCutoff))

    print('Read %d aligned reads of length %d. '
          'Found %d significant locations.' %
          (len(alignedReads), genomeLength, len(significantOffsets)))

    if not significantOffsets:
        print('Exiting due to finding no significant locations.')
        sys.exit(2)

    # agreementOffsets is a dict of dicts, where each value in the nested
    # dict is a tuple of sets. The first holds offsets that disagree, the
    # second those that agree. Both dicts are keyed by read id.
    agreementOffsets = defaultdict(
        lambda: defaultdict(
            lambda: (set(), set())
            )
        )

    DISAGREE, AGREE = (0, 1)

    # print('# counting agree/disagree at', ctime(time()), file=sys.stderr)
    for offset in significantOffsets:
        reads = readsAtOffset[offset]
        for read1 in reads:
            base1 = read1.base(offset)
            # assert base1
            for read2 in reads:
                if read2 is not read1:
                    base2 = read2.base(offset)
                    # assert base2
                    # The following line assumes AGREE is 1 & DISAGREE is 0.
                    agreement = int(base1 == base2)
                    agreementOffsets[read1][read2][agreement].add(offset)
                    agreementOffsets[read2][read1][agreement].add(offset)
    # print('# counted  agree/disagree at', ctime(time()), file=sys.stderr)

    readToNode = dict((read, Node(read)) for read in agreementOffsets)
    graphNodes = set()
    agreementCutoff = args.agreementCutoff

    # print('# making graph nodes at', ctime(time()), file=sys.stderr)
    for read1 in agreementOffsets:
        node1 = readToNode[read1]
        graphNodes.add(node1)
        for read2 in agreementOffsets[read1]:
            agreeCount = len(agreementOffsets[read1][read2][AGREE])
            disagreeCount = len(agreementOffsets[read1][read2][DISAGREE])
            total = agreeCount + disagreeCount
            assert total > 0
            if agreeCount / total >= agreementCutoff:
                node2 = readToNode[read2]
                node1.add(node2)
                node2.add(node1)
    # print('# made graph nodes at', ctime(time()), file=sys.stderr)

    if args.verbose:
        # Print the agree/disagree offsets and counts for all pairs of
        # reads.
        for read1 in agreementOffsets:
            print(read1)
            node1 = readToNode[read1]
            for read2 in agreementOffsets[read1]:
                node2 = readToNode[read2]
                print('  %s:\n    agree=%d (%s)\n    disagree=%d (%s)' % (
                    read2,
                    len(agreementOffsets[read1][read2][AGREE]),
                    sorted(agreementOffsets[read1][read2][AGREE]),
                    len(agreementOffsets[read1][read2][DISAGREE]),
                    sorted(agreementOffsets[read1][read2][DISAGREE])))
                if node2 in node1:
                    print('    Connected')

    # print('# computing ccs at', ctime(time()), file=sys.stderr)
    components = []
    for component in connectedComponents(graphNodes):
        offsets = componentOffsets(component)
        components.append((min(offsets), max(offsets), component))
    # print('# computed  ccs at', ctime(time()), file=sys.stderr)

    # Sort the components.
    def keyfunc(item):
        """
        Sort first on min offset then on component size.
        """
        return item[0], len(item[2])

    components.sort(key=keyfunc)

    if args.fastaDir:
        # Remove all pre-existing component-XXX.fasta and
        # component-XXX-aligned.fasta etc., files from the output
        # directory. This prevents us from doing a run that results in
        # (say) 6 component files and then doing a run that results in only
        # 5 and erroneously thinking that component-6.fasta is from the
        # second run.
        paths = list(map(str, chain(
            Path(args.fastaDir).glob('component-[0-9]*.fasta'),
            Path(args.fastaDir).glob('component-[0-9]*-aligned.fasta'),
            Path(args.fastaDir).glob('component-[0-9]*-base-frequencies.txt'),
            Path(args.fastaDir).glob('component-[0-9]*-consensus.fasta'),
            Path(args.fastaDir).glob('components-alignment.fasta'),
            Path(args.fastaDir).glob('components-alignment.txt'))))
        if paths:
            print('Removing %d pre-existing component .fasta/.txt file%s.' % (
                len(paths), '' if len(paths) == 1 else 's'))
            list(map(unlink, paths))

        componentCountWidth = int(log10(len(components))) + 1
        genomeLengthWidth = int(log10(genomeLength)) + 1
        readCountWidth = int(log10(len(alignedReads))) + 1
        consensuses = []

        for componentIndex, (minOffset, maxOffset, component) in enumerate(
                components):
            # The plain reads.
            filename = join(args.fastaDir,
                            'component-%0*d.fasta' % (
                                componentCountWidth, componentIndex + 1))
            with open(filename, 'w') as fp:
                for node in component:
                    print(node.read.read.toString('fasta'), end='', file=fp)

            # The aligned reads (i.e., with gaps in their sequences to
            # align them to the original genome).
            filename = join(
                args.fastaDir,
                'component-%0*d-aligned.fasta' % (
                    componentCountWidth, componentIndex + 1))
            with open(filename, 'w') as fp:
                if args.addConsensusToAlignment:
                    # Write the full consensus genome.
                    print(genome.toString('fasta'), end='', file=fp)
                for node in component:
                    preGaps = '-' * node.read.offset
                    postGaps = '-' * (
                        genomeLength - node.read.offset - len(node.read.read))
                    newRead = Read(
                        node.read.read.id,
                        preGaps + node.read.read.sequence + postGaps)
                    assert len(newRead) == genomeLength
                    print(newRead.toString('fasta'), end='', file=fp)

            # The base frequencies of component reads.
            filename = join(
                args.fastaDir,
                'component-%0*d-base-frequencies.txt' % (
                    componentCountWidth, componentIndex + 1))
            nucleotides = set('ACGT')
            consensus = []
            with open(filename, 'w') as fp:
                for offset in range(minOffset, maxOffset):
                    counts = Counter()
                    for node in component:
                        base = node.read.base(offset)
                        if base in nucleotides:
                            counts[base] += 1
                    print('Location %*d: base counts %s' % (
                        genomeLengthWidth, offset + 1,
                        baseCountsToStr(counts)), file=fp)
                    total = sum(counts.values())

                    # Figure out the most frequent bases.
                    maxCount = -1
                    bases = []
                    for base, baseCount in counts.items():
                        if baseCount > maxCount:
                            maxCount = baseCount
                            bases = [base]
                        elif baseCount == maxCount:
                            bases.append(base)
                    bases = ''.join(sorted(bases))
                    try:
                        consensus.append(BASES_TO_AMBIGUOUS[bases])
                    except KeyError:
                        if bases:
                            print(
                                'For component %d we have no ambiguous '
                                'nucleotide symbol for bases %r. Location %d '
                                '(min/max = %d/%d)' %
                                (componentIndex + 1, bases, offset + 1,
                                 minOffset + 1, maxOffset + 1),
                                file=sys.stderr)
                            # Put something into the consensus to serve as
                            # an alert.
                            consensus.append('[%s]' % (bases,))
                        else:
                            # No bases were found at this offset. This can
                            # happen if all reads in the component (and
                            # there may only be one) have (e.g.) an 'N' at
                            # this location.
                            consensus.append('N')

            # Write out the component consensus as FASTA.
            filename = join(
                args.fastaDir,
                'component-%0*d-consensus.fasta' % (
                    componentCountWidth, componentIndex + 1))
            consensusId = (
                'component-%0*d-consensus,reads=%d,locations=%d:%d,length=%d' %
                (componentCountWidth, componentIndex + 1, len(component),
                 minOffset + 1, maxOffset, maxOffset - minOffset))

            read = Read(consensusId, ''.join(consensus))
            with open(filename, 'w') as fp:
                print(read.toString('fasta'), end='', file=fp)

            consensuses.append((read, minOffset))

        # Write an alignment of the original consensus followed by all
        # component consensus sequences.
        filename = join(args.fastaDir, 'components-alignment.fasta')
        with open(filename, 'w') as fp:
            print(genome.toString('fasta'), end='', file=fp)
            for consensusRead, minOffset in consensuses:
                read = Read(
                    consensusRead.id,
                    ('-' * minOffset +
                     consensusRead.sequence +
                     '-' * (genomeLength - minOffset - len(consensusRead))))
                print(read.toString('fasta'), end='', file=fp)

        # See how well each component consensus fits with the original
        # (consensus) genome. That can be used to get a feel for whether
        # the components come from that genome or possibly some other.
        genomeSequence = genome.sequence
        filename = join(args.fastaDir, 'components-alignment.txt')
        with open(filename, 'w') as fp:
            for componentIndex, (consensusInfo, componentInfo) in enumerate(
                    zip(consensuses, components)):
                (consensus, minOffset) = consensusInfo
                consensusSequence = consensus.sequence
                (minOffset_, maxOffset, component) = componentInfo
                # Sanity check.
                assert minOffset == minOffset_

                agreeCount = 0
                for index in range(len(consensusSequence)):
                    if consensusSequence[index] in AMBIGUOUS[
                            genomeSequence[index + minOffset]]:
                        agreeCount += 1
                print('Component %*d (read count %*d, locations %*d-%*d, '
                      'length %*d): %*d of %*d (%.2f%%) bases match' %
                      (componentCountWidth, componentIndex + 1,
                       readCountWidth,
                       len(component),
                       genomeLengthWidth, minOffset + 1,
                       genomeLengthWidth, maxOffset,
                       genomeLengthWidth, maxOffset - minOffset,
                       genomeLengthWidth, agreeCount,
                       genomeLengthWidth, len(consensusSequence),
                       agreeCount / len(consensusSequence) * 100.0),
                      file=fp)

    # Write a textual summary of the components.
    for componentIndex, (minOffset, maxOffset, component) in enumerate(
            components):
        componentSize = len(component)
        print('Component %d (%d read%s) locations: %d-%d' % (
            componentIndex + 1, componentSize,
            '' if componentSize == 1 else 's', minOffset + 1, maxOffset))
        if args.verbose:
            for read in component:
                print('  ', read)

    # print('# making JSON at', ctime(time()), file=sys.stderr)
    if args.componentJSON:
        from json import dump
        result = []
        with open(args.componentJSON, 'w') as fp:
            for (minOffset, maxOffset, component) in components:
                componentReads = []
                for node in component:
                    alignedRead = node.read
                    componentReads.append(
                        (alignedRead.read.id, alignedRead.read.sequence,
                         alignedRead.offset))
                result.append(componentReads)
            dump(result, fp)
        del result
    # print('# made   JSON at', ctime(time()), file=sys.stderr)

    print('Found %d nodes in %d components. Agreement cutoff %.3f.' % (
        sum([len(component[-1]) for component in components]),
        len(components), agreementCutoff))

    if args.show:
        if args.additionalLocations:
            # These are 1-based.
            fields = args.additionalLocations.split()
            additionalOffsetsTitle = fields.pop(0).replace('_', ' ')
            additionalOffsets = list(
                map(lambda location: int(location) - 1, fields))
            assert all(offset >= 0 for offset in additionalOffsets)
        else:
            additionalOffsets = additionalOffsetsTitle = None

        title = (
            'Connected component composition<br>%d significant locations'
            % len(significantOffsets))
        plotComponents(
            components, args.readCategoryRegex, readCategoryRegexNames,
            genomeLength, significantOffsets, 'connected-components.html',
            title, additionalOffsets, additionalOffsetsTitle)

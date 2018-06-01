#!/usr/bin/env python

from __future__ import division, print_function

import sys
import plotly
import plotly.graph_objs as go
import re
import colorlover
from math import log10
from collections import defaultdict, Counter

from data.data import addCommandLineOptions, parseCommandLineOptions
from dark.dna import AMBIGUOUS
from dark.fasta import FastaReads

# Make a reverse version of AMBIGUOUS.
BASES_TO_AMBIGUOUS = {}
for symbol, bases in AMBIGUOUS.items():
    BASES_TO_AMBIGUOUS[''.join(sorted(bases))] = symbol


def connectedComponentsByOffset(significantReads):
    """
    Yield sets of reads that are connected according to what significant
    offsets they cover (the nucleotides at those offsets are irrelevant).

    @param significantReads: A C{set} of C{AlignedRead} instances, all of
        which cover at least one significant offset.
    """
    while significantReads:
        element = significantReads.pop()
        component = {element}
        offsets = set(element.significantOffsets)
        addedSomething = True
        while addedSomething:
            addedSomething = False
            these = set()
            for read in significantReads:
                if offsets.intersection(read.significantOffsets):
                    addedSomething = True
                    these.add(read)
                    offsets.update(read.significantOffsets)
            if these:
                significantReads.difference_update(these)
                component.update(these)
        yield component, offsets


def consistentReadSets(component, threshold):
    """
    Yield sets of reads that are consistent according to what nucleotides they
    have at their significant offsets.

    @param component: A C{set} of C{AlignedRead} instances.
    @param threshold: A C{float} indicating what fraction of read's nucleotides
        must be identical (to those already in the component) for it to be
        allowed to join a growing component.
    """
    def key(read):
        """
        We'll sort the avaialable reads by the number of significant offsets
        they have, then by start offset in the genome.
        """
        return (len(read.significantOffsets), read.offset)

    while component:
        reads = sorted(component, key=key, reverse=True)
        read0 = reads[0]
        these = {read0}
        nucleotides = defaultdict(Counter)
        for offset in read0.significantOffsets:
            nucleotides[offset][read0.base(offset)] += 1
        rejected = set()

        # Add all the reads that agree exactly at all the offsets in the
        # set so far.
        for read in reads[1:]:
            nucleotidesIfAccepted = []
            for offset in read.significantOffsets:
                base = read.base(offset)
                if offset in nucleotides:
                    if base not in nucleotides[offset]:
                        # Not an exact match. Reject this read for now.
                        rejected.add(read)
                        break
                nucleotidesIfAccepted.append((offset, base))
            else:
                # We didn't break, so add this read.
                for offset, base in nucleotidesIfAccepted:
                    nucleotides[offset][base] += 1
                component.remove(read)
                these.add(read)

        # Add in the first-round rejects that have a high enough threshold.
        # We do this before we add the bases of the rejects to nucleotides
        # because we don't want to pollute nucleotides with a bunch of
        # bases from first-round rejects (because their bases could
        # overwhelm the bases of the first round acceptees).
        acceptedRejects = set()
        for read in rejected:
            total = matching = 0
            for offset in read.significantOffsets:
                base = read.base(offset)
                if offset in nucleotides:
                    total += 1
                    if base in nucleotides[offset]:
                        matching += 1
            if total and matching / total >= threshold:
                # Looks good!
                acceptedRejects.add(read)
                component.remove(read)
                these.add(read)

        # Add the nucleotides of the second-round acceptances.
        for accepted in acceptedRejects:
            for offset in accepted.significantOffsets:
                nucleotides[offset][accepted.base(offset)] += 1

        component.difference_update(these)
        yield these, nucleotides


def nucleotidesToStr(nucleotides):
    """
    @param nucleotides: A C{defaultdict(Counter)} instance, keyed
        by offset, with nucleotides keying the Counters.
    """
    result = []
    for offset in sorted(nucleotides):
        result.append(
            '%d: %s' % (offset, baseCountsToStr(nucleotides[offset])))
    return ', '.join(result)


def baseCountsToStr(counts):
    """
    @param counts: A C{Counter} instance.
    """
    return ' '.join([
        ('%s:%d' % (base, counts[base])) for base in sorted(counts)])


def plotComponents(components, categoryRegexs, categoryNames, genomeLength,
                   significantOffsets, outFile, title, additionalOffsets,
                   additionalOffsetsTitle, show):
    """
    Show all connected components in a plot, as well as the significant
    locations and any additional locations.
    """
    offsetsForLevel = []
    categoryRegexs = categoryRegexs or []
    categoryNames = categoryNames or ['unnamed']
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
            id_ = node.read.id
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
    plotly.offline.plot(fig, filename=outFile, show_link=False, auto_open=show)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Find which reads agree and disagree with one another.')

    parser.add_argument(
        '--verbose', action='store_true', default=False,
        help='Print verbose textual output showing read connections.')

    parser.add_argument(
        '--reducedFilename',
        help=('The filename to write reduced (significant location only) '
              'reads to.'))

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
              'indicates to which genome a read belongs (when this is known, '
              'as is the case with simulated data.'))

    parser.add_argument(
        '--readCategoryRegexNames', nargs='*',
        help=('Specify names for the read categories defined in the regular '
              'expressions in --readCategoryRegex. If not specified, the '
              'regular expressions (i.e., their strings) will be used as the '
              'category names. If specified, must be the same length as '
              '--readCategoryRegex.'))

    parser.add_argument(
        '--agreementThreshold', type=float, default=0.9,
        help=('Only reads with agreeing nucleotides at at least this fraction '
              'of the significant sites they have in common will be '
              'considered connected.'))

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
              'being written for each component, due to --fastaDir being '
              'specified).'))

    addCommandLineOptions(parser)
    args = parser.parse_args()

    if args.readCategoryRegexNames:
        if not args.readCategoryRegex:
            print('--readCategoryRegex must be specified if '
                  '--readCategoryRegexNames is used.', file=sys.stderr)
            sys.exit(1)
        if len(args.readCategoryRegexNames) != len(args.readCategoryRegex):
            print('--readCategoryRegex and --readCategoryRegexNames must be '
                  'of equal length.', file=sys.stderr)
            sys.exit(1)

        readCategoryRegexNames = args.readCategoryRegexNames
    else:
        readCategoryRegexNames = list(args.readCategoryRegex or [])

    if args.runSpades and not args.fastaDir:
        print('--fastaDir must be specified if --runSpades is used.',
              file=sys.stderr)
        sys.exit(1)

    (genome, alignedReads, readCountAtOffset,
     baseCountAtOffset, readsAtOffset,
     significantOffsets) = parseCommandLineOptions(args, True)

    genomeLength = len(genome)

    print('Read %d aligned reads of length %d. '
          'Found %d significant locations.' %
          (len(alignedReads), genomeLength, len(significantOffsets)))

    if not significantOffsets:
        print('Exiting due to finding no significant locations.')
        sys.exit(2)

    if args.reducedFilename:
        allGaps = '-' * len(significantOffsets)
        def unwanted(read):
            return (None if read.sequence == allGaps else read)
        FastaReads(args.fastaFile).filter(
            keepSites=significantOffsets).filter(
                modifier=unwanted).save(args.reducedFilename)

    significantReads = set(read for read in alignedReads
                           if read.significantOffsets)

    componentCounts = []
    componentOffsets = []
    for i, (component, offsets) in enumerate(
            connectedComponentsByOffset(significantReads)):
        componentLength = len(component)
        print('component %d (%d reads)' % (i, componentLength))
        print('  offsets', sorted(offsets))
        for read in sorted(component):
            print('  ', read)

        subComponentCounts = []
        for j, (consistentReads, nucleotides) in enumerate(
                consistentReadSets(component, args.agreementThreshold)):
            subComponentCounts.append(len(consistentReads))
            print('    Consistent sub-graph %d (%d reads)' %
                  (j, len(consistentReads)))
            print('      nucleotides', nucleotidesToStr(nucleotides))
            for read in sorted(consistentReads):
                print('      ', read)

        assert componentLength == sum(subComponentCounts)
        subComponentCounts = sorted(subComponentCounts, reverse=True)
        print('  sub-component lengths', subComponentCounts)
        componentCounts.append(subComponentCounts)
        componentOffsets.append((min(offsets), max(offsets)))

    # Print a final component and sub-component count summary.
    print('\nFinal length summary of %d components' % len(componentCounts))
    totalReads = 0
    for i, (subComponentCounts, subComponentOffsets) in enumerate(
            zip(componentCounts, componentOffsets)):
        componentCount = sum(subComponentCounts)
        totalReads += componentCount
        print('component %d (%d reads) from %d to %d: %r' % (
            i, componentCount, subComponentOffsets[0], subComponentOffsets[1],
            subComponentCounts), end='')
        if len(subComponentCounts) > 1:
            print(' ratio %.2f' %
                  (subComponentCounts[0] / subComponentCounts[1]))
        else:
            print()

    print('%d reads (of %d) put into %d components.' % (
        totalReads, len(alignedReads), len(componentCounts)))


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

    # plotComponents(
    #     components, args.readCategoryRegex, readCategoryRegexNames,
    #     genomeLength, significantOffsets, 'connected-components.html',
    #     title, additionalOffsets, additionalOffsetsTitle, args.show)

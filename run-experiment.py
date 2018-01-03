#!/usr/bin/env python

from __future__ import division, print_function

import argparse
import sys
from time import time, ctime
from os import mkdir
from os.path import exists, join
from subprocess import run, PIPE, CalledProcessError
from random import choice
from json import load


class Executor(object):
    """
    Log and execute shell commands and store their stdout/stderr in an output
    directory.
    """
    def __init__(self, dryRun, logfp):
        self._dryRun = dryRun
        self._logFp = logfp

    def execute(self, command):
        """
        Use the shell to execute a command. Write stdout and stderr (if any)
        to the log file.
        """
        # Can't have newlines in a command given to the shell.
        command = command.replace('\n', ' ').strip()

        if self._dryRun:
            print('$ ' + command, file=sys.stderr)
            return

        start = time()
        print('# Start command at', ctime(start), file=self._logFp)
        print('$ ' + command, file=self._logFp)

        try:
            result = run(command, check=True, stdout=PIPE,
                         stderr=PIPE, shell=True, universal_newlines=True)
        except CalledProcessError as e:
            print('CalledProcessError:', e, file=self._logFp)
            print('STDOUT:\n%s' % e.stdout, file=self._logFp)
            print('STDERR:\n%s' % e.stderr, file=self._logFp)
            if e.stderr:
                print(e.stderr, file=sys.stderr)
            raise

        if result.stdout:
            print('STDOUT:\n%s' % result.stdout, file=self._logFp)

        if result.stderr:
            print('STDERR:\n%s' % result.stderr, file=self._logFp)

        stop = time()
        print('# Stop command at', ctime(stop), file=self._logFp)
        print('# Elapsed = %f seconds' % (stop - start), file=self._logFp)
        print(file=self._logFp)


def analyzeComponents(filename):
    with open(filename) as fp:
        components = load(fp)

    for component in components:
        size = len(component)
        genome1Count = 0
        for readId, readSequence, offset in component:
            if readId.startswith('genome-1-read-'):
                genome1Count += 1
        genome1Fraction = genome1Count / size
        print('Component split %.2f / %.2f (%d reads)' %
              (genome1Fraction, 1.0 - genome1Fraction, size))


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Run a dual infection experiment')

parser.add_argument(
    '--outputDir', required=True,
    help='The output directory. Must not already exist.')

genomeGroup = parser.add_mutually_exclusive_group(required=True)

genomeGroup.add_argument(
    '--genomeFile',
    help='The FASTA file containing the original genome to work from')

genomeGroup.add_argument(
    '--genomeLength', type=int, default=100,
    help='If a random initial genome is to be made, this will be its length')

parser.add_argument(
    '--verbose', action='store_true', default=False,
    help='Print verbose textual output.')

parser.add_argument(
    '--dryRun', action='store_true', default=False,
    help='If specified, simply print the actions that would be taken.')

parser.add_argument(
    '--genome2MutationRate', type=float, default=0.075,
    help='The per-base mutation rate to use to create genome 2')

parser.add_argument(
    '--readCount', type=int, default=100,
    help=('The number of reads to create for both genomees (only used '
          'if --genome-1-read-count or genome-2-read-count are not given).'))

parser.add_argument(
    '--minReads', type=int, default=5,
    help=('The minimum number of reads that must cover a location for it '
          'to be considered significant.'))

parser.add_argument(
    '--mixedReadsFile',
    help=('A FASTA file containing reads that are presumed to be from a '
          'multiple infection. If specified, reads from the original genome '
          'will not be created, a second genome will not be created, and '
          'reads from that second genome will also not be created. Processing '
          'will proceed directly to finding significant locations and then '
          'connected component analysis.'))

parser.add_argument(
    '--readCategoryRegex', nargs='*',
    help=('Specify regular expressions that match read ids for the '
          'purpose of plotting components and showing their catgeory '
          'composition. The category will normally be something that '
          'indicates to which genome a read belongs. This is only '
          '(potentially) useful when --mixedReadsFile is used and you have a '
          'way of distinguishing between the reads from the different '
          'infections. If --mixedReadsFile is not given, this argument is '
          'ignored.'))

parser.add_argument(
    '--readCategoryRegexNames', nargs='*',
    help=('Specify names for the read categories defined in the regular '
          'expressions in --readCategoryRegex. If not specified, the '
          'regular expressions (i.e., their strings) will be used as the '
          'category names. If specified, must be the same length as '
          '--readCategoryRegex. If --mixedReadsFile is not given, this '
          'argument is ignored.'))

parser.add_argument(
    '--homogeneousCutoff', type=float, default=0.85,
    help=('If the most common nucleotide at a location occurs more than '
          'this fraction of the time (i.e., amongst all reads that cover '
          'the location) then the locaion will be considered homogeneous '
          'and therefore uninteresting.'))

parser.add_argument(
    '--agreementCutoff', type=float, default=0.9,
    help=('Only reads that agree at least this much will be considered '
          'connected by connected-components.py.'))

parser.add_argument(
    '--genome1ReadCount', type=int,
    help='The number of reads to create for genome 1')

parser.add_argument(
    '--genome2ReadCount', type=int,
    help='The number of reads to create for genome 2')

parser.add_argument(
    '--minReadLength', type=int, default=10,
    help='The minimum length read to create')

parser.add_argument(
    '--maxReadLength', type=int, default=100,
    help=('The maximum length read to create. Defaults to the genome '
          'length'))

parser.add_argument(
    '--readMutationRate', type=float, default=0.0,
    help='The per-base mutation rate to use')

parser.add_argument(
    '--meanReadLength', type=float, default=100.0,
    help='The mean read length')

parser.add_argument(
    '--sdReadLength', type=float, default=10.0,
    help='The standard deviation of read length')

parser.add_argument(
    '--showComponents', action='store_true', default=False,
    help='If specified, show the connected component figure interactively.')

parser.add_argument(
    '--showBaseFrequencies', action='store_true', default=False,
    help='If specified, show the base frequencies figure interactively.')

parser.add_argument(
    '--additionalLocations',
    help=('Specify some additional locations to highlight. The argument '
          'must be a string name of the locations followed by space '
          'separated genome locations (1-based). Underscores in the '
          'name will be replaced by spaces.'))

args = parser.parse_args()
outputDir = args.outputDir

# We'll be using the shell to exec things and I don't want to have to quote
# all paths etc., so let's not allow whitespace in the output directory name.
# We should probably check for shell metacharacters too, but life is short.
if any(s.isspace() for s in outputDir):
    print('The output directory (%r) cannot have whitespace in its '
          'name. Exiting.' % outputDir, file=sys.stderr)
    sys.exit(1)

if exists(outputDir):
    if not args.dryRun:
        print('Output directory %r already exists. Exiting.' % outputDir,
              file=sys.stderr)
        sys.exit(2)
else:
    if not args.dryRun:
        mkdir(outputDir)

verbose = '--verbose' if args.verbose else ''
logfp = (sys.stderr if args.dryRun else
         open(join(outputDir, 'LOG'), 'w', buffering=1))

print('Command line arguments:\n  %s\n' % args, file=logfp)

genome1Fasta = join(outputDir, 'genome-1.fasta')
genome2Fasta = join(outputDir, 'genome-2.fasta')

genome1Alignment = join(outputDir, 'genome-1-alignment.fasta')
genome2Alignment = join(outputDir, 'genome-2-alignment.fasta')
mixedAlignment = join(outputDir, 'mixed-alignment.fasta')

executor = Executor(args.dryRun, logfp)

if args.genomeFile:
    executor.execute('cp %s %s' % (args.genomeFile, genome1Fasta))
else:
    if args.genomeLength < 1:
        print('Random initial genome length must be > 0.', file=sys.stderr)
        sys.exit(3)
    print('Writing random starting genome of length %d to %s' %
          (args.genomeLength, genome1Fasta), file=logfp)
    with open(genome1Fasta, 'w') as fp:
        print('>random-genome-1-length-%d\n%s' % (
            args.genomeLength,
            ''.join([choice('ACGT')
                     for _ in range(args.genomeLength)])), file=fp)

if args.mixedReadsFile:
    # If a mixed (i.e., presumed multiple infection) reads file has been
    # given, we skip various steps. Make sure incompatible options have not
    # been specified. Copy the mixed reads file to the expected location. (We
    # could instead just set mixedAlignment = args.mixedReadsFile.)
    executor.execute('cp %s %s' % (args.mixedReadsFile, mixedAlignment))
    if args.readCategoryRegex:
        readCategoryRegexArg = ' '.join(
            ["'%s'" % regex for regex in args.readCategoryRegex])
    else:
        # Match all reads in a single catch-all category.
        readCategoryRegexArg = '.'

    if args.readCategoryRegexNames:
        readCategoryRegexNamesArg = ' '.join(
            ["'%s'" % regex for regex in args.readCategoryRegexNames])
    else:
        # Give the reads a generic category name.
        readCategoryRegexNamesArg = "'Mixed infection'"
else:
    # Make a second genome using the given mutation rate.
    executor.execute('./mutate-reads.py --rate %s %s < %s > %s' %
                     (args.genome2MutationRate, verbose,
                      genome1Fasta, genome2Fasta))

    for info in [
            {
                'genomeAlignment': genome1Alignment,
                'genomeFasta': genome1Fasta,
                'genomeNumber': 1,
                'readCount': args.genome1ReadCount or args.readCount,
            },
            {
                'genomeAlignment': genome2Alignment,
                'genomeFasta': genome2Fasta,
                'genomeNumber': 2,
                'readCount': args.genome2ReadCount or args.readCount,
            }]:
        genomeNumber = info['genomeNumber']
        executor.execute(
            '''
            ./create-reads.py
                --maxReadLength %d
                --minReadLength %d
                --meanLength %d
                --sdLength %d
                --rate %s
                --idPrefix genome-%d-read-
                --count %d
                --alignReads
               < %s > %s
            ''' % (
                args.maxReadLength,
                args.minReadLength,
                args.meanReadLength,
                args.sdReadLength,
                args.readMutationRate,
                genomeNumber,
                info['readCount'],
                info['genomeFasta'],
                info['genomeAlignment']))

    executor.execute('cat %s %s > %s' %
                     (genome1Alignment, genome2Alignment, mixedAlignment))

    # We want 2 read id categories and we know what they are because we
    # just made the reads (see the --idPrefix arguments above in the calls
    # to create-reads.py). We need to take care with quoting here as these
    # will be part of a string that will be parsed by the shell.
    readCategoryRegexArg = "'^genome-1-' '^genome-2-'"
    readCategoryRegexNamesArg = "'Genome 1' 'Genome 2'"

connectedComponentJSON = join(outputDir, 'connected-components.json')
showComponents = '--show' if args.showComponents else ''
additionalLocations = ("--additionalLocations '%s'" % args.additionalLocations
                       if args.additionalLocations else '')

executor.execute('''
    ./connected-components.py
        --outFile %s
        --genomeFile %s
        --minReads %d
        --homogeneousCutoff %f
        --componentJSON %s
        --readCategoryRegex %s
        --readCategoryRegexNames %s
        --agreementCutoff %f
        --fastaDir %s
        --addConsensusToAlignment
        %s
        %s
        < %s > %s
''' % (
    join(outputDir, 'connected-components.html'),
    genome1Fasta,
    args.minReads,
    args.homogeneousCutoff,
    connectedComponentJSON,
    readCategoryRegexArg,
    readCategoryRegexNamesArg,
    args.agreementCutoff,
    outputDir,
    additionalLocations,
    showComponents,
    mixedAlignment,
    join(outputDir, 'connected-components.out')))

if not args.mixedReadsFile:
    # This analysis isn't aware of the read category regex list yet. It
    # only works when we made the mixed reads.
    analyzeComponents(connectedComponentJSON)

if args.showBaseFrequencies:
    executor.execute('''
        ./significant-base-frequencies.py
            --genomeFile %s
            --minReads %d
            --homogeneousCutoff %f
            --outFile %s
            --show
            < %s
    ''' % (
        genome1Fasta,
        args.minReads,
        args.homogeneousCutoff,
        join(outputDir, 'significant-base-frequencies.html'),
        mixedAlignment
    ))

logfp.close()

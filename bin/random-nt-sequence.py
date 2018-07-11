#!/usr/bin/env python

from random import choice
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Print a sequence of randomly chosen nucleotides to stdout.')

parser.add_argument(
    'n', type=int, default=75, help='The length of the sequence to print')

args = parser.parse_args()

print(''.join(choice('ACGT') for _ in range(args.n)))

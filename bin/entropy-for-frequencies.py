#!/usr/bin/env python

"""
Given category (label) frequencies on the command line, print the
entropy (base 2).

E.g., if there are 12 things in the first category, 4 in the second,
and 5 in the third, run

$ entropy.py 12 4 5
1.409975018822061
"""

import sys
from mid.entropy import entropy2

if len(sys.argv) < 2:
    print('Usage: %s freq1 [freq2 ...]' % sys.argv[0], file=sys.stderr)
    sys.exit(1)

try:
    counts = list(map(int, sys.argv[1:]))
except ValueError:
    print('%s: counts must be integers' % sys.argv[0], file=sys.stderr)
    sys.exit(2)

labels = []
for symbol, count in enumerate(counts):
    labels.extend([symbol] * count)

print(entropy2(labels))

def baseCountsToStr(counts):
    """
    Convert base counts to a string.

    @param counts: A C{Counter} instance.
    @return: A C{str} representation of nucleotide counts at an offset.
    """
    return ' '.join([
        ('%s:%d' % (base, counts[base])) for base in sorted(counts)])


def nucleotidesToStr(nucleotides, prefix=''):
    """
    Convert offsets and base counts to a string.

    @param nucleotides: A C{defaultdict(Counter)} instance, keyed
        by C{int} offset, with nucleotides keying the Counters.
    @param prefix: A C{str} to put at the start of each line.
    @return: A C{str} representation of the offsets and nucleotide
        counts for each.
    """
    result = []
    for offset in sorted(nucleotides):
        result.append(
            '%s%d: %s' % (prefix, offset,
                          baseCountsToStr(nucleotides[offset])))
    return '\n'.join(result)

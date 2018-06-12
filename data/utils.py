from contextlib import contextmanager
from pysam import AlignmentFile


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


def commonest(counts, drawFp=None, drawMessage=None):
    """
    Return the key of the Counter instance that is the most common.

    @param counts: A C{Counter} instance.
    @param drawFp: A file pointer to write information about draws (if any) to.
    @param drawMessage: A C{str} message to write to C{drawFp}. If the string
        contains '%(baseCounts)s' that will be replaced by a string
        representation of the base counts (in C{counts}) obtained from
        C{baseCountsToStr}. If not, the base count info will be printed after
        the message.
    """
    c = counts.most_common()

    # Detect draws & print a message if we find one.
    if len(c) > 1 and c[0][1] == c[1][1] and drawFp and drawMessage:
        bases = baseCountsToStr(counts)
        if drawMessage.find('%(baseCounts)s') > -1:
            print(drawMessage % {'baseCounts': bases}, file=drawFp)
        else:
            print('%s\n%s' % (drawMessage, bases), file=drawFp)

    return c[0][0]


@contextmanager
def samfile(filename):
    """
    A context manager to open and close a SAM/BAM file.

    @param filename: A C{str} file name to open.
    """
    f = AlignmentFile(filename)
    yield f
    f.close()

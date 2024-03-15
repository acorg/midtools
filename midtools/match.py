def _pct(a, b):
    """
    What percent of a is b?

    @param a: a numeric value.
    @param b: a numeric value.
    @return: the C{float} percentage.
    """
    return 100.0 * a / b if b else 0.0


def _pp(mesg, count, len1, len2=None):
    """
    Format a message followed by an integer count and a percentage (or
    two, if the sequence lengths are unequal).

    @param mesg: a C{str} message.
    @param count: a numeric value.
    @param len1: the C{int} length of sequence 1.
    @param len2: the C{int} length of sequence 2. If not given, will
        default to C{len1}.
    @return: A C{str} for printing.
    """
    if count == 0:
        return "%s: %d" % (mesg, count)
    else:
        len2 = len2 or len1
        if len1 == len2:
            return "%s: %d/%d (%.2f%%)" % (mesg, count, len1, _pct(count, len1))
        else:
            return (
                "%s: %d/%d (%.2f%%) of sequence 1, "
                "%d/%d (%.2f%%) of sequence 2)"
                % (mesg, count, len1, _pct(count, len1), count, len2, _pct(count, len2))
            )


def matchToString(dnaMatch, read1, read2, strict=False, indent="", offsets=None):
    """
    Format a DNA match as a string.

    @param dnaMatch: A C{dict} returned by C{compareDNAReads}.
    @param strict: If C{True}, ambiguous nucleotide symbols were not allowed
        to match.
    @return: A C{str} describing the match.
    """
    match = dnaMatch["match"]
    identicalMatchCount = match["identicalMatchCount"]
    ambiguousMatchCount = match["ambiguousMatchCount"]
    gapMismatchCount = match["gapMismatchCount"]
    gapGapMismatchCount = match["gapGapMismatchCount"]
    nonGapMismatchCount = match["nonGapMismatchCount"]

    if offsets:
        len1 = len2 = len(offsets)
    else:
        len1, len2 = map(len, (read1, read2))

    result = []
    append = result.append

    append(_pp("%sExact matches" % indent, identicalMatchCount, len1, len2))
    append(_pp("%sAmbiguous matches" % indent, ambiguousMatchCount, len1, len2))
    if ambiguousMatchCount and identicalMatchCount:
        anyMatchCount = identicalMatchCount + ambiguousMatchCount
        append(_pp("%sExact or ambiguous matches" % indent, anyMatchCount, len1, len2))
    mismatchCount = gapMismatchCount + gapGapMismatchCount + nonGapMismatchCount
    append(_pp("%sMismatches" % indent, mismatchCount, len1, len2))
    conflicts = "conflicts or ambiguities" if strict else "conflicts"
    append(
        _pp(
            "%s  Not involving gaps (i.e., %s)" % (indent, conflicts),
            nonGapMismatchCount,
            len1,
            len2,
        )
    )
    append(
        _pp(
            "%s  Involving a gap in one sequence" % indent, gapMismatchCount, len1, len2
        )
    )
    append(
        _pp(
            "%s  Involving a gap in both sequences" % indent,
            gapGapMismatchCount,
            len1,
            len2,
        )
    )

    for read, key in zip((read1, read2), ("read1", "read2")):
        append("%s  Id: %s" % (indent, read.id))
        length = len(read)
        append("%s    Length: %d" % (indent, length))
        gapCount = len(dnaMatch[key]["gapOffsets"])
        append(_pp("%s    Gaps" % indent, gapCount, length))
        if gapCount:
            append(
                "%s    Gap locations (1-based): %s"
                % (
                    indent,
                    ", ".join(
                        map(
                            lambda offset: str(offset + 1),
                            sorted(dnaMatch[key]["gapOffsets"]),
                        )
                    ),
                )
            )
        ambiguousCount = len(dnaMatch[key]["ambiguousOffsets"])
        append(_pp("%s    Ambiguous" % indent, ambiguousCount, length))
        extraCount = dnaMatch[key]["extraCount"]
        if extraCount:
            append(_pp("%s    Extra nucleotides at end" % indent, extraCount, length))

    return "\n".join(result)

#!/usr/bin/env python

from __future__ import division, print_function

import plotly
import plotly.graph_objs as go
from json import load


def sampleIdKey(id_):
    """
    Get a sort key for a sample id.

    @param id_: The C{str} id of a sample. This will normally take the form
        of DA222 (letters followed by some digits), but may also have trailing
        non-digits, as in VK211-G.
    @raise ValueError: If id_ contains a digit that is not followed by
        characters that are all also digits.
    @return: A 3-C{tuple} containing the initial C{str} part of id_ then
        the C{int} suffix, if any (else C{0}), then the trailing C{str}, if
        any.
    """
    # This is stolen from pyhbv/samples.py
    prefix = ""
    value = 0
    trailing = ""

    state = "PREFIX"

    for offset, c in enumerate(id_):
        if state == "PREFIX":
            if c.isdigit():
                state = "VALUE"
                value = 10 * value + int(c)
            else:
                prefix += c
        elif state == "VALUE":
            if c.isdigit():
                value = 10 * value + int(c)
            else:
                state = "TRAILING"
                trailing += c
        else:
            trailing += c

    return prefix, value, trailing


def plotScores(sampleData, outfile):
    """
    Plot the sorted base scores for the samples in C{sampleData}.
    """
    data = []
    maxLength = -1

    for sampleName in sorted(sampleData, key=sampleIdKey):
        values = sampleData[sampleName]["values"]
        length = len(values)
        if length > maxLength:
            maxLength = length
        initialZeroCount = 0
        for value in values:
            if value == 0:
                initialZeroCount += 1
            else:
                break
        data.append(
            go.Scatter(
                x=list(range(1 + initialZeroCount, length + 1)),
                y=values[initialZeroCount:],
                # mode='lines+markers',
                mode="markers",
                name=sampleName,
                text=sampleData[sampleName]["text"][initialZeroCount:],
                marker=dict(
                    size=4,
                    # color=color,
                    line=dict(
                        width=0,
                        # color=color,
                    ),
                ),
            )
        )

    xmargin = max(1, int(maxLength * 0.01))
    xaxis = {
        "title": "Rank",
        "range": (-xmargin, maxLength + xmargin),
    }

    yaxis = {
        "title": "Significant location nucleotide frequency entropy",
    }

    layout = go.Layout(xaxis=xaxis, yaxis=yaxis)
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename=outfile, show_link=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Show multiple sorted significant genome location "
            "nucleotide frequencies for a set of aligned reads."
        ),
    )

    parser.add_argument(
        "json",
        nargs="+",
        metavar="file.json",
        help=(
            "The JSON files containing the sorted location values. These "
            "are produced using the --valuesFile option when running "
            "significant-base-frequencies.py"
        ),
    )

    parser.add_argument(
        "--outfile",
        default="multiple-significant-base-frequencies.html",
        help="The filename to store the resulting HTML.",
    )

    args = parser.parse_args()

    data = {}
    for filename in args.json:
        with open(filename, "r") as fp:
            x = load(fp)
            data[x["sampleName"]] = {
                "values": x["values"],
                "text": x["text"],
            }

    plotScores(data, args.outfile)

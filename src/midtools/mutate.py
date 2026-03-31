from random import uniform, choice

from dark.reads import Read


def mutateRead(read: Read, rate: float, alphabet: str = "ACGT") -> list[int]:
    """
    Mutate the bases of a read.

    @param read: A C{dark.reads.Read} instance.
    @param rate: A C{float} per-base mutation rate.
    @param alphabet: An iterable of single character strings. If characters are
        repeated this will increase their probability of being selected to
        replace existing characters in the read.
    @return: A C{list} of the offsets where C{read} was mutated.
    """
    newSequence = []
    mutatedOffsets = []
    for offset, base in enumerate(read.sequence):
        if uniform(0.0, 1.0) <= rate:
            # Ignore the possibility of an infinite loop due to Rosencrantz
            # & Guildenstern Are Dead class bad luck or to 'alphabet' only
            # containing one thing (the base we are trying to mutate away
            # from).
            while True:
                newBase = choice(alphabet)
                if newBase != base:
                    break
            newSequence.append(newBase)
            mutatedOffsets.append(offset)
        else:
            newSequence.append(base)

    read.sequence = "".join(newSequence)

    return mutatedOffsets

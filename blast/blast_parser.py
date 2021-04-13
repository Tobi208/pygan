from typing import List, Tuple, TextIO

Read = List[Tuple[str, float]]
Accessions = List[str]


def parse(file: str, blast_format: str = 'tab', top_score_percent: float = 1) -> List[Accessions]:
    """
    Read lines of file, extract reads containing accessions and bit scores.
    Filter accessions in each read by the top score percentage.
    Assumes that reads are continuous and read ids don't need to be stored.

    :param file: filepath
    :param blast_format: blast format (tab, xml, pairwise)
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :return: list of accessions per read filtered by top score percentage
    """

    with open(file, 'r') as f:

        accessions = None
        if blast_format == 'tab':
            accessions = parse_tab(f, top_score_percent)
        elif blast_format == 'xml':
            pass
        elif blast_format == 'pairwise':
            pass

    return accessions


def parse_tab(f: TextIO, top_score_percent: float) -> List[Accessions]:
    """
    Read lines of file in tab format, extract reads containing accessions and bit scores.
    Filter accessions in each read by the top score percentage.
    Assumes that reads are continuous and read ids don't need to be stored.

    qseqid = 0
    sseqid = 1
    bitscore = 2

    :param f: file in tab format
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :return: list of accessions per read filtered by top score percentage
    """

    all_accessions = []
    read = None
    read_id = None
    line = f.readline()

    while line:
        line = line.strip('\n').split('\t')
        next_id = line[0]
        # first read
        if not read:
            read_id = next_id
            read = []
        # new read
        if next_id != read_id:
            # flush last read
            all_accessions.append(filter_by_top_score(read, top_score_percent))
            read = []
            read_id = next_id
        # expand read
        accession = line[1][:-2]
        bit_score = float(line[2])
        read.append((accession, bit_score))
        # iterate
        line = f.readline()

    # flush last read
    all_accessions.append(filter_by_top_score(read, top_score_percent))

    return all_accessions


def filter_by_top_score(read: Read, top_score_percent: float) -> Accessions:
    """
    Filter accessions in read by the top score percentage.
    An accession within the percentage of the top score stays in the read.

    Example: top score = 50, top_score_percent = 0.1: 47 remains, 43 is discarded.

    :param read: read to be filtered
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :return list of accessions with scores >= top_score_percent of top score
    """

    top_score = 0
    for _, s in read:
        if s > top_score:
            top_score = s

    bound = top_score - top_score * top_score_percent
    return [accession for accession, score in read if not score < bound]

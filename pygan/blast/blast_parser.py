from typing import List, Tuple, Dict


def parse_filter(file: str, top_score_percent: float, tab_map: Dict[str, int]) -> Tuple[List[List[str]], List[str]]:
    """
    Read lines of file in tab format, extract reads containing accessions and bit scores.
    Filter accessions in each read by the top score percentage.
    Assumes that reads are continuous and read ids don't need to be stored.

    :param file: filepath
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :param tab_map: contains a mapping of which column qseqid, sseqid and bitscore are in
    :return: list of accessions per read filtered by top score percentage, list of read ids
    """

    qseqid = tab_map['qseqid']
    sseqid = tab_map['sseqid']
    bitscore = tab_map['bitscore']

    reads = []
    read_ids = []
    read = None
    read_id = None

    with open(file, 'r') as f:
        line = f.readline()
        while line:
            line = line.strip('\n').split('\t')
            next_id = line[qseqid]
            # first read
            if not read:
                read_id = next_id
                read = []
            # new read
            if next_id != read_id:
                # flush last read
                reads.append(filter_by_top_score(read, top_score_percent))
                read_ids.append(read_id)
                read = []
                read_id = next_id
            # expand read
            accession = line[sseqid][:-2]
            bit_score = float(line[bitscore])
            read.append((accession, bit_score))
            # iterate
            line = f.readline()
        # flush last read
        reads.append(filter_by_top_score(read, top_score_percent))
        read_ids.append(read_id)

    return reads, read_ids


def filter_by_top_score(read: List[Tuple[str, float]], top_score_percent: float) -> List[str]:
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

from typing import List
from pygan.tree.phylo_tree import PhyloTree


# megan rank ids -> scientific names map
RANKS = {
    '0': 'unspecified',
    '1': 'kingdom',
    '2': 'phylum',
    '3': 'class',
    '4': 'order',
    '5': 'family',
    '90': 'varietas',
    '98': 'genus',
    '99': 'species group',
    '100': 'species',
    '101': 'subspecies',
    '127': 'domain'
}

TAXON_KEY = 0
NAME_KEY = 1
RANK_KEY = 3


def map_names(file: str, tree: PhyloTree):
    """
    Read lines of file, extract tax_id to name mapping and apply it to tree

    :param file: filepath
    :param tree: phylo tree to be filled with names
    """
    parse(readlines(file), tree)


def readlines(file: str) -> List[List[str]]:
    """
    Read lines of a file and return list of first two tab separated entries per line

    :param file: filepath
    :return: list of first two tab separated entries per line
    """
    with open(file, 'r') as f:
        return [line.split('\t') for line in f.read().splitlines()]


def parse(lines: List[List[str]], tree: PhyloTree):
    """
    Fill tree nodes with names and ranks for corresponding tax_ids

    :param lines: list of [tax_id, name, ..., rank, ...]
    :param tree: phylo tree to be filled with names and ranks
    """
    nodes = tree.nodes
    for line in lines:
        taxon_id = int(line[TAXON_KEY])
        if taxon_id in nodes:
            node = nodes[taxon_id]
            node.name = line[NAME_KEY]
            rank_id = line[RANK_KEY]
            if rank_id in RANKS:
                node.rank = RANKS[rank_id]
            else:
                node.rank = RANKS['0']

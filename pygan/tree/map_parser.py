from typing import List
from pygan.tree.phylo_tree import PhyloTree


# megan rank ids -> scientific names map
ranks = {
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
    taxon = 0
    name = 1
    rank = 3

    for line in lines:
        if line[taxon] in nodes:
            node = nodes[line[taxon]]
            node.name = line[name]
            rank_id = line[rank]
            if rank_id in ranks:
                node.rank = ranks[rank_id]
            else:
                node.rank = ranks['0']

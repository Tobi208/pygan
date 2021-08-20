from typing import Dict, List, Optional


class PhyloNode:
    """
    Primitive phylogenetic tree node

    Contains a tax_id, name, rank, path from root, pointer to its parent, list of children and
    and reads indicates the number of reads mapped to this node or is a list of their ids.
    """

    def __init__(self):
        self.tax_id: Optional[int] = None
        self.name: Optional[str] = None
        self.name_with_rank: Optional[str] = None
        self.rank: Optional[str] = None
        self.path: Optional[str] = None
        self.path_with_rank: Optional[str] = None
        self.reads: List[str] = []
        self.parent: Optional[PhyloNode] = None
        self.children: List[PhyloNode] = []

    def to_string(self, show_path: bool, list_reads: bool, show_rank: bool):
        """
        Print a phylogenetic node with various options

        :param show_path: prefix path from root to node
        :param list_reads: display read ids mapped to node (else display number of reads)
        :param show_rank: add an abbreviation of the rank to the name
        :return:
        """
        reads = ','.join(self.reads) if list_reads else str(len(self.reads))
        spaced = '\t' + reads + '\n'
        if show_path and show_rank:
            return self.path_with_rank + spaced
        elif show_path and not show_rank:
            return self.path + spaced
        elif show_rank:
            return self.name_with_rank + spaced
        else:
            return self.name + spaced


class PhyloTree:
    """
    Primitive phylogenetic tree

    Contains only pointer to its root and tax_id to node map
    """

    _RANK_ABBREV: Dict[Optional[str], str] = {
        None: '',
        'unspecified': '',
        'kingdom': 'k',
        'phylum': 'p',
        'class': 'c',
        'order': 'o',
        'family': 'f',
        'varietas': 'v',
        'genus': 'g',
        'species group': 'sg',
        'species': 's',
        'subspecies': 'ss',
        'domain': 'd'
    }

    def __init__(self):
        self.root: Optional[PhyloNode] = None
        self.nodes: Dict[int, PhyloNode] = {}

    def clear_reads(self):
        """
        Remove all mapped reads from the tree
        """
        for node in self.nodes.values():
            node.reads.clear()

    def completed_mapping(self):
        """
        After mapping is completed, parse additional info to string in nodes
        """
        self._add_rank_generate_path(self.root, '', '')

    def _add_rank_generate_path(self, node: PhyloNode, path: str, path_with_rank: str):
        """
        Add rank as a prefix and path to current node

        :param node: current node
        :param path: accumulated path to current node
        :param path_with_rank: accumulated path to current node with ranks prefixed
        :return:
        """
        # attempt to prefix rank
        if self._RANK_ABBREV[node.rank]:
            node.name_with_rank = self._RANK_ABBREV[node.rank] + '__' + node.name
        else:
            node.name_with_rank = node.name
        # accumulate path
        path += '/' + node.name
        path_with_rank += '/' + node.name_with_rank
        node.path = path
        node.path_with_rank = path_with_rank
        # propagate
        for child in node.children:
            self._add_rank_generate_path(child, path, path_with_rank)

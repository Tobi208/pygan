from typing import Dict, List, Optional, Union


rank_abbrev = {
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


class PhyloTree:
    """
    Primitve phylogenetic tree

    Contains only pointer to its root and tax_id to node map and some flags
    """

    def __init__(self):
        self.root: Optional[PhyloNode] = None
        self.nodes: Dict[int, PhyloNode] = {}
        self.reads_type: type = list
        self.ranks_prefixed: bool = False
        self.paths_generated: bool = False

    def convert_to_num_reads(self):
        for node in self.nodes.values():
            node.reads = len(node.reads)
        self.reads_type = int

    def add_rank_abbrev_to_name(self):
        for node in self.nodes.values():
            if rank_abbrev[node.rank]:
                node.name = rank_abbrev[node.rank] + '__' + node.name
        self.ranks_prefixed = True
        if self.paths_generated:
            self.generate_paths()

    def generate_paths(self):
        def paths_dfs(node: PhyloNode, path):
            path += '/' + node.name
            node.path = path
            for child in node.children:
                paths_dfs(child, path)
        paths_dfs(self.root, '')
        self.paths_generated = True


class PhyloNode:
    """
    Primitive phylogentic tree node

    Contains a tax_id, name, rank, path from root, pointer to its parent, list of children and
    and reads indicates the number of reads mapped to this node or is a list of their ids.
    """

    def __init__(self):
        self.tax_id: Optional[int] = None
        self.name: Optional[str] = None
        self.rank: Optional[str] = None
        self.path: Optional[str] = None
        self.reads: Union[List[str], int] = []
        self.parent: Optional[PhyloNode] = None
        self.children: List[PhyloNode] = []

    def to_string(self, show_path, list_reads):
        s = ''
        if show_path:
            s += self.path
        else:
            s += self.name
        s += '\t'
        if list_reads:
            s += ','.join(self.reads)
        else:
            s += str(self.reads)
        s += '\n'
        return s
from typing import Dict, List, Optional


class PhyloTree:
    """
    Primitve phylogenetic tree

    Contains only pointer to its root and tax_id to node map
    """

    def __init__(self):
        self.root: Optional[PhyloNode] = None
        self.nodes: Dict[int, PhyloNode] = {}


class PhyloNode:
    """
    Primitive phylogentic tree node

    Contains only a tax_id, name, rank, pointer to its parent and list of children.
    Reads indicates the number of reads mapped to this node.
    """

    def __init__(self):
        self.tax_id: Optional[int] = None
        self.name: Optional[str] = None
        self.rank: Optional[str] = None
        self.reads: List[str] = []
        self.parent: Optional[PhyloNode] = None
        self.children: List[PhyloNode] = []

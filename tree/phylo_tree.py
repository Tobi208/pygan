class PhyloTree:
    """
    Primitve phylogenetic tree

    Contains only pointer to its root and tax_id to node map
    """

    def __init__(self):
        self.root = None
        self.nodes = {}


class PhyloNode:
    """
    Primitive phylogentic tree node

    Contains only a tax_id, name, rank, pointer to its parent and list of children.
    Reads indicates the number of reads mapped to this node.
    """

    def __init__(self):
        self.tax_id = None
        self.name = None
        self.rank = None
        self.reads = None
        self.parent = None
        self.children = []

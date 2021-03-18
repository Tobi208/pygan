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

    Contains only a tax_id, name, pointer to its parent and list of children
    """

    def __init__(self):
        self.tax_id = None
        self.name = None
        self.parent = None
        self.children = []

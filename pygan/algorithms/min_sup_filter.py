from pygan.tree.phylo_tree import PhyloTree, PhyloNode


def apply(tree: PhyloTree, min_support: int):
    """
    Applies the minimum support filter to a phylogenetic tree.
    Every node with a reads count lower than the specified limit donates its reads to its parent.

    :param tree: phylogenetic tree with #reads mapped to nodes
    :param min_support: minimum support limit nodes are required to satisfy
    """
    root = tree.root
    if root:
        min_sup_dfs(root, min_support)


def min_sup_dfs(node: PhyloNode, min_support: int):
    """
    Enforces a minimum support limit bottom up on a phylogenetic tree by traversing the tree depth first.

    :param node: current node in the phylogenetic tree (call first on root)
    :param min_support: minimum support limit nodes are required to satisfy
    """
    for child in node.children:
        min_sup_dfs(child, min_support)
    # if a node's reads count is below the min sup limit
    # push its reads upwards
    if 0 < len(node.reads) < min_support and node.parent:
        node.parent.reads[:] += node.reads
        node.reads.clear()
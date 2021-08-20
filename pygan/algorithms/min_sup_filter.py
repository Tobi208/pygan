from pygan.tree.phylo_tree import PhyloTree, PhyloNode

from typing import List

ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species group', 'species', 'subspecies',
         'varietas', 'domain']
major_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def apply(tree: PhyloTree, min_support: int, exclude: List[str], only_major: bool = False):
    """
    Applies the minimum support filter to a phylogenetic tree.
    Every node with a reads count lower than the specified limit donates its reads to its parent.

    :param tree: phylogenetic tree with #reads mapped to nodes
    :param only_major: only major ranks are allowed to retain reads
    :param min_support: minimum support limit nodes are required to satisfy
    :param exclude: ranks that the filter should not be applied to
    """
    root = tree.root
    if root:
        if only_major:
            min_sup_dfs_to_major(root, min_support, exclude)
        else:
            min_sup_dfs(root, min_support, exclude)


def min_sup_dfs_to_major(node: PhyloNode, min_support: int, exclude: List[str]):
    """
    Enforces a minimum support limit bottom up on a phylogenetic tree by traversing the tree depth first.
    Only nodes with majors ranks are allowed to retain reads.

    :param node: current node in the phylogenetic tree (call first on root)
    :param min_support: minimum support limit nodes are required to satisfy
    :param exclude: ranks that the filter should not be applied to
    """
    for child in node.children:
        min_sup_dfs_to_major(child, min_support, exclude)
        # if a node's reads count is below the min sup limit or its not a major rank
        # push its reads upwards
    if (node.rank not in major_ranks or 0 < len(node.reads) < min_support) and node.parent and node.rank not in exclude:
        node.parent.reads[:] += node.reads
        node.reads.clear()


def min_sup_dfs(node: PhyloNode, min_support: int, exclude: List[str]):
    """
    Enforces a minimum support limit bottom up on a phylogenetic tree by traversing the tree depth first.

    :param node: current node in the phylogenetic tree (call first on root)
    :param min_support: minimum support limit nodes are required to satisfy
    :param exclude: ranks that the filter should not be applied to
    """
    for child in node.children:
        min_sup_dfs(child, min_support, exclude)
    # if a node's reads count is below the min sup limit
    # push its reads upwards
    if 0 < len(node.reads) < min_support and node.parent and node.rank not in exclude:
        node.parent.reads[:] += node.reads
        node.reads.clear()

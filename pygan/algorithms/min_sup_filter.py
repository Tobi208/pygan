from math import ceil
from typing import Dict
from pygan.tree.phylo_tree import PhyloTree, PhyloNode


ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species group', 'species', 'subspecies',
         'varietas', 'domain']
major_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


def apply(tree: PhyloTree, min_support: int, only_major: bool = False):
    """
    Applies the minimum support filter to a phylogenetic tree.
    Every node with a reads count lower than the specified limit donates its reads to its parent.

    :param tree: phylogenetic tree with #reads mapped to nodes
    :param only_major: only major ranks are allowed to retain reads
    :param min_support: minimum support limit nodes are required to satisfy
    """
    root = tree.root
    if root:
        if only_major:
            min_sup_dfs_to_major(root, min_support)
        else:
            min_sup_dfs(root, min_support)


def min_sup_dfs_to_major(node: PhyloNode, min_support: int):
    """
    Enforces a minimum support limit bottom up on a phylogenetic tree by traversing the tree depth first.
    Only nodes with majors ranks are allowed to retain reads.

    :param node: current node in the phylogenetic tree (call first on root)
    :param min_support: minimum support limit nodes are required to satisfy
    """
    for child in node.children:
        min_sup_dfs_to_major(child, min_support)
        # if a node's reads count is below the min sup limit or its not a major rank
        # push its reads upwards
    if (node.rank not in major_ranks or 0 < len(node.reads) < min_support) and node.parent:
        node.parent.reads[:] += node.reads
        node.reads.clear()


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


def project_to_rank(tree: PhyloTree, rank: str):
    """
    Project all reads to nodes of a certain rank.
    Since reads are projected, specific read ids are discarded
    and number of reads is used instead.
    Readsare projected downwards by percentage.

    :param tree: phylogenetic tree with reads mapped to nodes
    :param rank: name of rank to project reads to
    """
    tree.convert_to_num_reads()
    root = tree.root
    if root:
        project_dfs_up(root, rank, {})
        project_dfs_down(root, rank)


def project_dfs_up(node: PhyloNode, rank: str, push_up: Dict[int, bool]):
    """
    Push reads upwards until specified rank is hit.
    Then sum reads upwards to enable downward projection.

    :param node: current node in the phylogenetic tree
    :param rank: rank to project reads to
    :param push_up: dictionary mapping whether to push up reads of a node (by tax_id)
    """
    # check if leaf is specified rank
    if not node.children:
        push_up[node.tax_id] = not is_rank(node, rank)
    # traveserse dfs if not a leaf
    else:
        for child in node.children:
            project_dfs_up(child, rank, push_up)
    # if not root
    parent = node.parent
    if parent:
        # propagate reads upwards
        parent.reads += node.reads
        # check if parent be deleting reads
        # once rank is met, stop deleting reads
        if parent.tax_id in push_up:
            push_up[parent.tax_id] = push_up[node.tax_id] and push_up[parent.tax_id]
        else:
            push_up[parent.tax_id] = push_up[node.tax_id] and not is_rank(parent, rank)
        # if node should delete reads, delete them
        if push_up[node.tax_id]:
            node.reads = 0


def project_dfs_down(node: PhyloNode, rank: str):
    """
    Push reads downwards until specified rank is hit.
    Distribute reads by shares of child reads.

    :param node: current node in the phylogenetic tree
    :param rank: rank to project reads to
    """
    # gather amount of reads in children
    total_child_reads = sum(child.reads for child in node.children)
    # distribute node reads by shares of child reads
    # and traverse dfs
    for child in [child for child in node.children if child.reads > 0]:
        share = child.reads / total_child_reads
        child.reads = ceil(share * node.reads)
        project_dfs_down(child, rank)
    # delete nodes if above specified rank
    if not is_rank(node, rank):
        node.reads = 0


def is_rank(node: PhyloNode, rank: str):
    """

    :param node: node in phylogenetic tree
    :param rank: rank to compare to
    :return: if node rank matches rank
    """
    return node.rank != 'unspecified' and node.rank == rank

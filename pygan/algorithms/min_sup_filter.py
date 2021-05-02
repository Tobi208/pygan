from pygan.tree.phylo_tree import PhyloTree, PhyloNode

major_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

ranks = {
    'kingdom': 1,
    'phylum': 2,
    'class': 3,
    'order': 4,
    'family': 5,
    'varietas': 90,
    'genus': 98,
    'species group': 99,
    'species': 100,
    'subspecies': 101,
    'domain': 127,
}


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


def project_to_rank(tree: PhyloTree, rank_name: str):
    """
    Project all reads to nodes of a certain rank.
    Since reads are projected, specific read ids are discarded
    and number of reads is used instead.
    Reads above the specified rank are projected downwards by percentage.

    :param tree: phylogenetic tree with #reads mapped to nodes
    :param rank_name: name of rank to project reads to
    """
    tree.convert_to_num_reads()
    root = tree.root
    if root:
        project_dfs_up(root, ranks[rank_name])
        project_dfs_down(root, ranks[rank_name])


def project_dfs_up(node: PhyloNode, rank: int):
    """
    Push reads upwards until specified rank is hit.
    Then sum reads upwards to enable downward projection.

    :param node: current node in the phylogenetic tree
    :param rank: rank to project reads to
    """
    for child in node.children:
        project_dfs_up(child, rank)
    if node.parent:
        node.parent.reads += node.reads
    if node.rank == 'unspecified' or ranks[node.rank] > rank:
        node.reads = 0


def project_dfs_down(node: PhyloNode, rank: int):
    """

    :param node: current node in the phylogenetic tree
    :param rank: rank to project reads to
    """
    if node.rank == 'unspecified' or ranks[node.rank] < rank:
        total_child_reads = sum(child.reads for child in node.children)
        for child in node.children:
            if total_child_reads > 0:
                share = child.reads / total_child_reads
                child.reads = share * node.reads
            project_dfs_down(child, rank)
        node.reads = 0

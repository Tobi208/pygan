from typing import List, Optional, Dict
from math import ceil

from pygan.tree.phylo_tree import PhyloTree, PhyloNode


def project_proportional(tree: PhyloTree, rank: str):
    """
    Collect reads from nodes below the target rank,
    then project reads from nodes above the target rank proportionally downwards.

    :param tree: phylo tree to project reads along
    :param rank: rank to project reads to
    """

    proportional_up(tree.root, rank, False)
    proportional_down(tree.root)


def proportional_up(node: PhyloNode, rank: str, below: bool):
    """
    Collect reads from nodes below the target rank
    and mark nodes that can have their reads projected downwards later.

    :param node: current node
    :param rank: target rank to project reads to
    :param below: flag if current node is below the target rank
    """

    # children will modify marked and sum of (projectable) reads
    node.marked = False
    node.sum_of_reads = 0

    # propagate
    for child in node.children:
        proportional_up(child, rank, below or node.rank == rank)

    # if below rank, pass reads upwards
    if below:
        node.parent.reads += node.reads
        node.reads.clear()
    # if target rank, begin upwards marking if reads are present
    elif node.rank == rank and len(node.reads) > 0:
        node.marked = True

    # if branch is projectable, propagate upwards
    if node.marked and node.parent is not None:
        node.sum_of_reads += len(node.reads)
        node.parent.sum_of_reads += node.sum_of_reads
        node.parent.marked = True


def proportional_down(node: PhyloNode):
    """
    Project reads of node proportionally downwards if possible.

    :param node: current node
    """

    # only project to marked children and if there is anything to project
    prj_children = [child for child in node.children if child.marked and child.sum_of_reads > 0]
    total_children_reads = sum(child.sum_of_reads for child in prj_children)

    # pass reads proportionally downwards
    if len(node.reads) > 0 and total_children_reads > 0:
        i = 0
        for child in prj_children:
            share = ceil(child.sum_of_reads / total_children_reads * len(node.reads))
            j = i + int(share)
            child.reads += node.reads[i:j]
            i = j
        if i < len(node.reads) - 1:
            prj_children[-1].reads += node.reads[i:]
        node.reads.clear()

    # propagate
    for child in prj_children:
        proportional_down(child)


def project_accession(tree: PhyloTree, rank: str, reads: List[List[int]], read_ids: List[str], cluster_degree: int):
    """
    Project reads to a target rank by remapping them to lower level taxons according to their accessions
    or by pushing them upwards if they are below the target rank.

    :param tree: phylo tree
    :param rank: target rank for projection
    :param reads: list of potential taxons for each read
    :param read_ids: list of read ids
    :param cluster_degree: degree of clustering of low level taxons
    """

    accession_up(tree.root, rank, None)
    read_map = {read_id: taxids for read_id, taxids in zip(read_ids, reads)}
    accession_down(tree.nodes, tree.root, rank, read_map, cluster_degree)


def accession_up(node: PhyloNode, rank: str, target: Optional[PhyloNode]):
    """
    Collect reads from nodes below the target rank.

    :param node: current node
    :param rank: target rank to project reads to
    :param target: node that children below rank should push their reads to
    """

    # if projection target exists, project to target
    if target is not None:
        target.reads += node.reads
        node.reads.clear()
    # if node is of target rank, set as target
    elif node.rank == rank:
        target = node

    # propagate
    for child in node.children:
        accession_up(child, rank, target)


def accession_down(nodes: Dict[int, PhyloNode], node: PhyloNode, rank: str,
                   read_map: Dict[str, List[int]], cluster_degree: int):
    """
    Project reads to a target rank by remapping them to lower level taxons according to their accessions.

    :param nodes: dict of nodes in the phylo tree
    :param node: current node to have its reads projected
    :param rank: target rank of projection
    :param read_map: map of read ids to the read's mapped accessions
    :param cluster_degree: degree of clustering of low level taxons
    """

    # terminate recursion at target rank
    if node.rank == rank:
        return

    # reads that can not be projected and are retained by the node
    retain = []

    # attempt to project reads
    for read_id in node.reads:

        # determine hits
        taxids = read_map[read_id]
        hits = get_hits(nodes, rank, taxids)
        # if there are no suitable candidates, disregard
        if len(hits) == 0:
            retain.append(read_id)
            continue

        # cluster and evaluate hits
        clusters = get_clusters(nodes, hits, cluster_degree)
        best_hit = get_best_hit(clusters, hits)
        # reassign read
        nodes[best_hit].reads.append(read_id)

    # retain non-projectable reads
    node.reads = retain

    # propagate
    for child in node.children:
        accession_down(nodes, child, rank, read_map, cluster_degree)


def get_hits(nodes: Dict[int, PhyloNode], rank: str, taxids: List[int]) -> Dict[int, int]:
    """
    Collect taxons in their ancestors of the target rank.

    :param nodes: dict of nodes in the phylo tree
    :param rank: target rank of projection
    :param taxids: ids of low level taxons corresponding to a read
    :return: ancestors with the target rank and how many taxons map to them
    """

    hits = {}
    for taxid in taxids:
        if taxid not in nodes:
            continue
        hit = get_ancestor_of_rank(nodes[taxid], rank)
        # pigeonhole ancestors of taxons
        if not hit:
            continue
        if hit in hits:
            hits[hit] += 1
        else:
            hits[hit] = 1
    return hits


def get_ancestor_of_rank(node: PhyloNode, rank: str) -> Optional[int]:
    """
    Get ancestor of a node with a specified rank.

    :param node: current node
    :param rank: target rank of projection
    :return: ancestor of node with target rank
    """

    while node:
        if node.rank == rank:
            return node.tax_id
        node = node.parent
    return None


def get_clusters(nodes: Dict[int, PhyloNode], hits: Dict[int, int], cluster_degree: int) \
        -> Dict[PhyloNode, List[int]]:
    """
    Cluster hits to a specified degree.

    :param nodes: dict of nodes in the phylo tree
    :param hits: potential nodes to map the read to
    :param cluster_degree: degree to cluster the hits
    :return: clusters of hits
    """

    clusters = {}
    # pigeonhole ancestors of a specific degree of hits
    for taxid in hits.keys():
        cluster = nodes[taxid]
        for _ in range(cluster_degree):
            if cluster.parent is not None:
                cluster = cluster.parent
        if cluster in clusters:
            clusters[cluster].append(taxid)
        else:
            clusters[cluster] = [taxid]
    return clusters


def get_best_hit(clusters: Dict[PhyloNode, List[int]], hits: Dict[int, int]) -> int:
    """
    Determine the highest scoring hit of the highest scoring cluster

    :param clusters: clusters of hits
    :param hits: potential nodes to map the read to
    :return: best hit
    """

    # determine best cluster
    best_cluster = -1
    best_val = -1
    for cluster, cluster_hits in clusters.items():
        val = sum(hits[hit] for hit in cluster_hits)
        if val > best_val:
            best_val = val
            best_cluster = cluster
    # determine best hit in a cluster
    best_hit = -1
    best_val = -1
    for taxid in clusters[best_cluster]:
        if hits[taxid] > best_val:
            best_val = hits[taxid]
            best_hit = taxid
    return best_hit


def project_accession_proportional(tree: PhyloTree, rank: str,
                                   reads: List[List[int]], read_ids: List[str], cluster_degree: int):
    """
    Project reads to a target rank by remapping them to lower level taxons according to their accessions
    or by pushing them upwards if they are below the target rank. Then project remaining reads proportionally.

    :param tree: phylo tree
    :param rank: target rank for projection
    :param reads: list of potential taxons for each read
    :param read_ids: list of read ids
    :param cluster_degree: degree of clustering of low level taxons
    """

    project_accession(tree, rank, reads, read_ids, cluster_degree)
    project_proportional(tree, rank)

from typing import Dict, Tuple, List, Any
from pickle import dump, load
from time import time
from pygan.tree.phylo_tree import PhyloTree
from pygan.tree.newick_parser import get_phylo_tree
from pygan.tree.map_parser import map_names
from pygan.blast.blast_parser import parse_filter, parse_with_score, filter_by_top_score
from pygan.database.megan_map import get_accessions2taxonids
from pygan.algorithms.lca import compute_addresses, get_common_prefix
from pygan.algorithms.min_sup_filter import apply
from pygan.algorithms.projection import project_proportional, project_accession, project_accession_proportional


def run(tre_file: str, map_file: str, megan_map_file: str, blast_file: str,
        blast_map: Dict[str, int], top_score_percent: float, db_segment_size: int, db_key: str,
        ignore_ancestors: bool, min_support: int, only_major: bool, exclude: List[str],
        project_mode: str, project_rank: str, cluster_degree: int,
        out_file: str, prefix_rank: bool, show_path: bool, list_reads: bool):
    """
    Performs an LCA analysis

    LCA analysis takes user input from a blast file,
    maps its accessions to taxon ids via the Megan Map Database,
    determines the Lowest Common Ancestor for each read and
    maps the number of reads to a phylogenetic tree.
    Read projection can be applied afterward.
    Print the resulting taxonomy to a text file.

    :param tre_file: path to file containing phylogenetic tree
    :param map_file: path to file containing mapping of taxonomy id to scientific name and rank
    :param megan_map_file: path to file containing megan_map.db
    :param blast_file: path to file containing blast data
    :param blast_map: contains a mapping of which column qseqid, sseqid and bitscore are in
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :param db_segment_size: number of reads whose accessions are to be mapped via the database in chunks
    :param db_key: specific key to map accessions to (Taxonomy for NCBI, gtdb for GTDB)
    :param ignore_ancestors: flag whether to ignore ancestors in the LCA algorithm
    :param min_support: limit for the minimum support filter algorithm
    :param only_major: only major ranks are allowed to retain reads
    :param exclude: ranks that the minimum support filter should not be applied to
    :param project_mode: method of read projection to a target rank ('accession', 'proportional', 'mixed', '')
    :param project_rank: target rank to project reads to
    :param cluster_degree: degree to cluster hits by during the accession-based projection
    :param out_file: path to output file of results
    :param prefix_rank: add an abbreviation of the rank to the name
    :param show_path: show paths from root to nodes
    :param list_reads: display read ids mapped to nodes (else display number of reads)
    """

    print('starting lca analysis')
    lca_start = time()
    tree = parse_tree(tre_file, map_file)
    id2address, address2id = compute_lca_addresses(tree)
    reads, read_ids = parse_blast_filter(blast_file, top_score_percent, blast_map)
    mapped_reads = map_accessions(reads, megan_map_file, db_segment_size, db_key)
    map_lcas(tree, id2address, address2id, mapped_reads, read_ids, ignore_ancestors)
    project_reads_to_rank(project_mode, tree, project_rank, mapped_reads, read_ids, cluster_degree)
    apply_min_sup_filter(tree, min_support, exclude, only_major)
    write_results(tree, out_file, prefix_rank, show_path, list_reads)
    print('completed lca analysis in ' + timer(lca_start))


def parse_tree(tre_file: str, map_file: str) -> PhyloTree:
    """
    Parse a phylogenetic tree from a newick format file
    and map names and ranks from a map file to it

    :param tre_file: path to file containing phylogenetic tree
    :param map_file: path to file containing mapping of taxonomy id to scientific name and rank
    :return: phylogenetic tree with taxonomy ids, scientific names and ranks
    """
    t = time()
    tree = get_phylo_tree(tre_file)
    print('parsed tree in ' + timer(t))
    t = time()
    map_names(map_file, tree)
    tree.completed_mapping()
    print('mapped names and ranks in ' + timer(t))
    return tree


def compute_lca_addresses(tree: PhyloTree) -> Tuple[Dict, Dict]:
    """
    Compute a mapping of taxonomy id to its address in the tree and vice versa

    :param tree: phylogenetic tree
    :return: mapping of taxonomy id to its address in the tree and vice versa
    """
    t = time()
    id2address = {}
    address2id = {}
    compute_addresses(tree, id2address, address2id)
    print('computed addresses in ' + timer(t))
    return id2address, address2id


def parse_blast_filter(blast_file: str, top_score_percent: float, blast_map: Dict[str, int]) \
        -> Tuple[List[List[str]], List[str]]:
    """
    Parse a blast tab file and filter the accessions by top score

    :param blast_file: path to file containing blast data
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :param blast_map: contains a mapping of which column qseqid, sseqid and bitscore are in
    :return: list of accessions per read filtered by top score, list of read ids
    """
    t = time()
    reads_n_read_ids = parse_filter(blast_file, top_score_percent, blast_map)
    print('parsed blast in ' + timer(t))
    return reads_n_read_ids


def parse_blast_with_score(blast_file: str, blast_map: Dict[str, int]) \
        -> Tuple[List[List[Tuple[str, float]]], List[str]]:
    """
    Parse a blast tab file with scores

    :param blast_file: path to file containing blast data
    :param blast_map: contains a mapping of which column qseqid, sseqid and bitscore are in
    :return: list of accessions with score per read, list of read ids
    """
    t = time()
    reads_ws_n_read_ids = parse_with_score(blast_file, blast_map)
    print('parsed blast with score in ' + timer(t))
    return reads_ws_n_read_ids


def filter_reads_by_top_score(reads_ws: List[List[Tuple[Any, float]]], top_score_percent: float) -> List[List[Any]]:
    """
    Filter items in read by the top score percentage.
    An item within the percentage of the top score stays in the read.

    Example: top score = 50, top_score_percent = 0.1: 47 remains, 43 is discarded.

    :param reads_ws: list of reads with scores
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :return list of items with scores >= top_score_percent of top score
    """
    t = time()
    reads = [filter_by_top_score(read, top_score_percent) for read in reads_ws]
    print('filtered reads by top score in ' + timer(t))
    return reads


def map_accessions(reads: List[List[str]], megan_map_file: str, db_segment_size: int, db_key: str) -> List[List[int]]:
    """
    Retrieve taxonomy ids for every read from the Megan Map Database

    :param reads: list of accessions per read
    :param megan_map_file: path to file containing megan_map.db
    :param db_segment_size: number of reads whose accessions are to be mapped via the database in chunks
    :param db_key: specific key to map accessions to (Taxonomy for NCBI, gtdb for GTDB)
    :return: list of taxonomy ids per read
    """
    t = time()
    mapped_reads = []
    # compute indices for segmentation
    segments = [*range(0, len(reads), db_segment_size), len(reads)]
    for i in range(1, len(segments)):
        # group reads into chunks
        grouped_reads = reads[segments[i - 1]:segments[i]]
        # collect all accessions from a chunk of reads
        flattened_reads = [acc for read in grouped_reads for acc in read]
        # map accessions to taxons
        acc2id = get_accessions2taxonids(megan_map_file, flattened_reads, db_key)
        # dechunk reads again
        for read in grouped_reads:
            mapped_reads.append([acc2id[acc] for acc in read if acc in acc2id])
    print('mapped #reads: ' + str(len(reads)) + ' in ' + timer(t))
    return mapped_reads


def map_accessions_with_scores(reads_ws: List[List[Tuple[str, float]]],
                               megan_map_file: str, db_segment_size: int, db_key: str) \
        -> List[List[Tuple[int, float]]]:
    """
    Retrieve taxonomy ids for every read from the Megan Map Database

    :param reads_ws: list of accessions with scores per read
    :param megan_map_file: path to file containing megan_map.db
    :param db_segment_size: number of reads whose accessions are to be mapped via the database in chunks
    :param db_key: specific key to map accessions to (Taxonomy for NCBI, gtdb for GTDB)
    :return: list of taxonomy ids with scores per read
    """
    t = time()
    mapped_reads_ws = []
    # compute indices for segmentation
    segments = [*range(0, len(reads_ws), db_segment_size), len(reads_ws)]
    for i in range(1, len(segments)):
        # group reads into chunks
        grouped_reads_ws = reads_ws[segments[i - 1]:segments[i]]
        # collect all accessions from a chunk of reads
        flattened_reads = [acc for read_ws in grouped_reads_ws for acc, _ in read_ws]
        # map accessions to taxons
        acc2id = get_accessions2taxonids(megan_map_file, flattened_reads, db_key)
        # dechunk reads again
        for read_ws in grouped_reads_ws:
            mapped_reads_ws.append([(acc2id[acc], score) for acc, score in read_ws if acc in acc2id])
    print('mapped #reads: ' + str(len(reads_ws)) + ' in ' + timer(t))
    return mapped_reads_ws


def map_lcas(tree: PhyloTree, id2address: Dict, address2id: Dict,
             reads: List[List[int]], read_ids: List[str], ignore_ancestors: bool):
    """
    Computes Lowest Common Ancestors for each read and maps it to the corresponding node in the phylogenetic tree

    :param tree: phylogenetic tree
    :param id2address: mapping of taxonomy id to its address in the tree
    :param address2id: mapping of a tree address to its taxonomy id
    :param reads: list of taxonomy ids per read
    :param read_ids: list of read ids corresponding to reads
    :param ignore_ancestors: use longest address or shortest address as reference
    """
    t = time()
    nodes = tree.nodes
    # map each read to a taxon
    for i, read in enumerate(reads):
        # by computing the common prefix on its mapped accessions
        common_prefix = get_common_prefix([
            id2address[taxonid] for taxonid in read
            if taxonid in id2address
        ], ignore_ancestors)
        # map read
        nodes[address2id[common_prefix]].reads.append(read_ids[i])
    print('computed LCAs in ' + timer(t))


def apply_min_sup_filter(tree: PhyloTree, min_support: int, exclude: List[str], only_major: bool):
    """
    Applies the minimum support filter to a phylogenetic tree

    :param tree: phylogenetic tree
    :param min_support: minimum support limit nodes are required to satisfy
    :param only_major: only major ranks are allowed to retain reads
    :param exclude: ranks that the filter should not be applied to
    """
    t = time()
    if min_support < 2 and not only_major:
        return
    apply(tree, min_support, exclude, only_major)
    print('applied min support filter in ' + timer(t))


def project_reads_to_rank(mode: str, tree: PhyloTree, rank: str,
                          mapped_reads: List[List[int]], read_ids: List[str], cluster_degree: int):
    """
    Attempt to project reads of nodes to only nodes with the target rank.
    Reads mapped below the target rank are pushed upwards.
    For reads mapped above the target rank, multiple heuristics are available.

    :param mode: projection mode: proportional, accession or mixed
    :param tree: phylo tree
    :param rank: target rank for projection
    :param mapped_reads: list of potential taxons for each read
    :param read_ids: list of read ids
    :param cluster_degree: degree of clustering of low level taxons
    """
    t = time()
    if mode == 'proportional':
        project_proportional(tree, rank)
    elif mode == 'accession':
        project_accession(tree, rank, mapped_reads, read_ids, cluster_degree)
    elif mode == 'mixed':
        project_accession_proportional(tree, rank, mapped_reads, read_ids, cluster_degree)
    else:
        return
    print('projected reads to rank in ' + timer(t))


def write_results(tree: PhyloTree, out_file: str, prefix_rank: bool, show_path: bool, list_reads: bool):
    """
    Writes the results of the lca analysis to a file

    :param tree: phylogenetic tree
    :param out_file: path to output file
    :param prefix_rank: add an abbreviation of the rank to the name
    :param show_path: show paths from root to nodes
    :param list_reads: display read ids mapped to nodes (else display number of reads)
    """
    t = time()
    result = [node.to_string(show_path, list_reads, prefix_rank) for node in tree.nodes.values() if node.reads]
    with open(out_file, 'w') as f:
        f.writelines(result)
    print('exported result in ' + timer(t))


def save_to_bin(obj: Any, file: str):
    """
    Save object to binary file with pickle

    :param obj: object to save to binary file
    :param file: output file
    """
    t = time()
    with open(file, 'wb') as f:
        dump(obj, f)
    print('saved to binary in ' + timer(t))


def load_from_bin(file: str) -> Any:
    """
    Load object from binary file with pickle

    :param file: binary file containing object
    :return: object from binary file
    """
    t = time()
    with open(file, 'rb') as f:
        obj = load(f)
    print('loaded from binary in ' + timer(t))
    return obj


def timer(t: float) -> str:
    """
    :param t: reference time
    :return: formatted time difference
    """
    return str(round(time() - t, 2))

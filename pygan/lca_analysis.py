from time import time
from typing import Dict, Tuple, List
from pygan.tree.phylo_tree import PhyloTree
from pygan.tree.newick_parser import get_phylo_tree
from pygan.tree.map_parser import map_names
from pygan.blast.blast_parser import parse_filter
from pygan.database.megan_map import get_accessions2taxonids
from pygan.algorithms.lca import compute_addresses, get_common_prefix
from pygan.algorithms.min_sup_filter import apply


def run(tre_file: str, map_file: str, megan_map_file: str, blast_file: str,
        blast_map: Dict[str, int], top_score_percent: float, db_segment_size: int, db_key: str,
        ignore_ancestors: bool, min_support: int, out_file: str):
    """
    Conducts an LCA analysis

    LCA analysis takes user input from a blast file,
    maps its accessions to taxon ids via the Megan Map Database,
    determines the Lowest Common Ancestor for each read and
    maps the number of reads to a phylogenetic tree.

    :param tre_file: path to file containing phyolgenetic tree
    :param map_file: path to file containing mapping of taxonomy id to scientific name and rank
    :param megan_map_file: path to file containing megan_map.db
    :param blast_file: path to file containing blast data
    :param blast_map: contains a mapping of which column qseqid, sseqid and bitscore are in
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :param db_segment_size: number of reads whoose accessions are to be mapped via the database in chunks
    :param db_key: specific key to map accessions to (Taxonomy for NCBI, gtdb for GTDB)
    :param ignore_ancestors: flag whether to ignore ancestors in the LCA algorithm
    :param min_support: limit for the minimum support filter algorithm
    :param out_file: path to output file of results
    """

    print('starting lca analysis')
    lca_start = time()
    tree = parse_tree(tre_file, map_file)
    id2address, address2id = compute_lca_addresses(tree)
    reads, read_ids = parse_blast_filter(blast_file, top_score_percent, blast_map)
    mapped_reads = map_accessions(reads, megan_map_file, db_segment_size, db_key)
    map_lcas(tree, id2address, address2id, mapped_reads, read_ids, ignore_ancestors)
    apply_min_sup_filter(tree, min_support)
    write_results(tree, out_file)
    print('completed lca analysis in ' + timer(lca_start))


def parse_tree(tre_file: str, map_file: str) -> PhyloTree:
    """
    Parse a pyhlogenetic tree from a newick format file
    and map names and ranks from a map file to it

    :param tre_file: path to file containing phyolgenetic tree
    :param map_file: path to file containing mapping of taxonomy id to scientific name and rank
    :return: phylogenetic tree with taxonomy ids, scientific names and ranks
    """
    t = time()
    tree = get_phylo_tree(tre_file)
    print('parsed tree in ' + timer(t))
    t = time()
    map_names(map_file, tree)
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


def parse_blast_filter(blast_file: str, top_score_percent: float, blast_map: Dict[str, int])\
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


def map_accessions(reads: List[List[str]], megan_map_file: str, db_segment_size: int, db_key: str) -> List[Tuple[int]]:
    """
    Retrieve taxonomy ids for every read from the Megan Map Database

    :param reads: list of accessions per read
    :param megan_map_file: path to file containing megan_map.db
    :param db_segment_size: number of reads whoose accessions are to be mapped via the database in chunks
    :param db_key: specific key to map accessions to (Taxonomy for NCBI, gtdb for GTDB)
    :return: list of taxonomy ids per read
    """
    t = time()
    mapped_reads = []
    segments = [*range(0, len(reads), db_segment_size), len(reads)]
    for i in range(1, len(segments)):
        grouped_reads = reads[segments[i - 1]:segments[i]]
        flattened_reads = [acc for read in grouped_reads for acc in read]
        acc2id = get_accessions2taxonids(megan_map_file, flattened_reads, db_key)
        for read in grouped_reads:
            mapped_reads.append(tuple(acc2id[acc] for acc in read if acc in acc2id))
    print('mapped #reads: ' + str(len(reads)) + ' in ' + timer(t))
    return mapped_reads


def map_lcas(tree: PhyloTree, id2address: Dict, address2id: Dict,
             reads: List[Tuple[int]], read_ids: List[str], ignore_ancestors: bool):
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
    for i, read in enumerate(reads):
        common_prefix = get_common_prefix([
            id2address[taxonid] for taxonid in read
            if taxonid in id2address
        ], ignore_ancestors)
        nodes[address2id[common_prefix]].reads.append(read_ids[i])
    print('computed LCAs in ' + timer(t))


def apply_min_sup_filter(tree: PhyloTree, min_support: int):
    """
    Applies the minimum support filter to a phylogenetic tree

    :param tree: phylogenetic tree
    :param min_support: minimum support limit nodes are required to satisfy
    """
    t = time()
    apply(tree, min_support)
    print('applied min support filter in ' + timer(t))


def write_results(tree: PhyloTree, out_file: str):
    """
    Writes the results of the lca analysis to a file

    :param tree: phylogenetic tree
    :param out_file: path to output file
    """
    t = time()
    result = '\n'.join([node.name + '\t' + str(len(node.reads)) for node in tree.nodes.values() if node.reads])
    with open(out_file, 'w') as f:
        f.write(result)
    print('exported result in: ' + timer(t))


def timer(t: float) -> str:
    """
    :param t: reference time
    :return: formatted time difference
    """
    return str(round(time() - t, 2))

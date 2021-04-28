from time import time
from pygan.tree.newick_parser import get_phylo_tree
from pygan.tree.map_parser import map_names
from pygan.blast.blast_parser import parse as blast_parse
from pygan.database.megan_map import get_accessions2taxonids
from pygan.algorithms.lca import compute_addresses, get_common_prefix
from pygan.algorithms.min_sup_filter import apply_min_sup_filter


def execute(tre_file: str, map_file: str, megan_map_file: str, blast_file: str,
            blast_format: str, top_score_percent: float, db_segment_size: int, db_key: str,
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
    :param blast_format: specific blast format (tab, xml, pairwise)
    :param top_score_percent: percentage in [0, 1] to filter accessions by
    :param db_segment_size: number of reads whoose accessions are to be mapped via the database in chunks
    :param db_key: specific key to map accessions to (Taxonomy for NCBI, gtdb for GTDB)
    :param ignore_ancestors: flag whether to ignore ancestors in the LCA algorithm
    :param min_support: limit for the minimum support filter algorithm
    :param out_file: path to output file of results
    """

    print('starting lca analysis')
    lca_start = time()

    t = time()
    tree = get_phylo_tree(tre_file)
    print('parsed tree in ' + timer(t))

    t = time()
    id2address = {}
    address2id = {}
    compute_addresses(tree, id2address, address2id)
    print('computed addresses in ' + timer(t))

    t = time()
    map_names(map_file, tree)
    print('mapped names and ranks in ' + timer(t))

    t = time()
    all_accessions = blast_parse(blast_file, blast_format, top_score_percent)
    print('parsed blast in ' + timer(t))

    t = time()
    segments = [*range(0, len(all_accessions), db_segment_size), len(all_accessions)]
    all_taxonids = []
    for i in range(1, len(segments)):
        grouped_accs = all_accessions[segments[i - 1]:segments[i]]
        flattened_accs = [acc for accs in grouped_accs for acc in accs]
        acc2id = get_accessions2taxonids(megan_map_file, flattened_accs, db_key)
        all_taxonids += [tuple(acc2id[acc] for acc in accs if acc in acc2id) for accs in grouped_accs]
    print('mapped #reads: ' + str(len(all_taxonids)) + ' in ' + timer(t))

    t = time()
    nodes = tree.nodes
    for taxonids in all_taxonids:
        common_prefix = get_common_prefix([
            id2address[taxonid] for taxonid in taxonids
            if taxonid in id2address
        ], ignore_ancestors)
        nodes[address2id[common_prefix]].reads += 1
    print('computed LCAs in ' + timer(t))

    t = time()
    apply_min_sup_filter(tree, min_support)
    print('applied min support filter in ' + timer(t))

    t = time()
    result = '\n'.join([node.name + '\t' + str(node.reads) for node in tree.nodes.values() if node.reads != 0])
    with open(out_file, 'w') as f:
        f.write(result)
    print('exported result in: ' + timer(t))

    print('completed lca analysis in ' + timer(lca_start))


def timer(t):
    return str(round(time() - t, 2))

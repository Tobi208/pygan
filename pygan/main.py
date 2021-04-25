from time import time
from pygan.tree.newick_parser import get_phylo_tree
from pygan.tree.map_parser import map_names
from pygan.blast.blast_parser import parse as blast_parse
from pygan.database.megan_map import fetch_all_taxonids
from pygan.algorithms.lca import compute_addresses, get_common_prefix
from pygan.algorithms.min_sup_filter import apply_min_sup_filter


def lca_analysis(tre_file: str, map_file: str, megan_map_file: str, blast_file: str,
                 blast_format: str, top_score_percent: float, ignore_ancestors: bool, min_support: int,
                 out_file: str):

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
    all_taxonids = fetch_all_taxonids(megan_map_file, all_accessions)
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


if __name__ == '__main__':
    lca_analysis(tre_file='../resources/ncbi.tre',
                 map_file='../resources/ncbi.map',
                 megan_map_file='../resources/megan-map-Jan2021.db',
                 blast_file='../resources/Alice01-1mio-Jan-2021.txt',
                 blast_format='tab', top_score_percent=0.1, ignore_ancestors=False, min_support=100,
                 out_file='../lca_analysis.txt')
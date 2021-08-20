from pygan.lca_analysis import *


def run_full_lca():
    """
    Example on how to run an LCA analysis with a single call.
    """
    run(tre_file='resources/ncbi.tre',
        map_file='resources/ncbi.map',
        megan_map_file='resources/megan-map-Jan2021.db',
        blast_file='resources/Alice01-1mio-Jan-2021.txt',
        blast_map={'qseqid': 0, 'sseqid': 1, 'bitscore': 2},
        top_score_percent=0.1,
        db_segment_size=10000, db_key='Taxonomy',
        ignore_ancestors=False, min_support=100, only_major=False,
        exclude=[], project_mode='', project_rank='', cluster_degree=0,
        out_file='lca_analysis.txt',
        prefix_rank=True,
        show_path=False,
        list_reads=False)


def script_lca():
    """
    Example on how to interactively script an LCA analysis.
    """
    tre_file = 'resources/ncbi.tre'
    map_file = 'resources/ncbi.map'
    megan_map_file = 'resources/megan-map-Jan2021.db'
    blast_file = 'resources/Alice01-1mio-Jan-2021.txt'
    blast_map = {'qseqid': 0, 'sseqid': 1, 'bitscore': 2}
    top_score_percent = 0.1
    db_segment_size = 10000
    db_key = 'Taxonomy'
    ignore_ancestors = False
    min_support = 50
    only_major = False
    exclude = ['genus']
    project_mode = 'proportional'
    project_rank = 'genus'
    cluster_degree = 2
    out_file = 'lca_analysis.txt'
    prefix_rank = True
    show_path = False
    list_reads = False

    tree = parse_tree(tre_file, map_file)
    # tree = load_from_bin('tree_post_lca.bin')
    id2address, address2id = compute_lca_addresses(tree)
    reads, read_ids = parse_blast_filter(blast_file, top_score_percent, blast_map)
    mapped_reads = map_accessions(reads, megan_map_file, db_segment_size, db_key)
    map_lcas(tree, id2address, address2id, mapped_reads, read_ids, ignore_ancestors)
    # save_to_bin(tree, 'tree_post_lca.bin')
    project_reads_to_rank(project_mode, tree, project_rank, mapped_reads, read_ids, cluster_degree)
    apply_min_sup_filter(tree, min_support, exclude, only_major)
    write_results(tree, out_file, prefix_rank, show_path, list_reads)


if __name__ == '__main__':
    run_full_lca()

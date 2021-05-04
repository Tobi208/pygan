from pygan.lca_analysis import *


def run_full_lca():
    run(tre_file='resources/ncbi.tre',
        map_file='resources/ncbi.map',
        megan_map_file='resources/megan-map-Jan2021.db',
        blast_file='resources/Alice01-1mio-Jan-2021.txt',
        blast_map={'qseqid': 0, 'sseqid': 1, 'bitscore': 2},
        top_score_percent=0.1,
        db_segment_size=10000, db_key='Taxonomy',
        ignore_ancestors=False, min_support=100, only_major=False,
        out_file='lca_analysis.txt',
        prefix_rank=True,
        show_path=False,
        list_reads=False)


def msf_only_major():
    run(tre_file='resources/ncbi.tre',
        map_file='resources/ncbi.map',
        megan_map_file='resources/megan-map-Jan2021.db',
        blast_file='resources/Alice01-1mio-Jan-2021.txt',
        blast_map={'qseqid': 0, 'sseqid': 1, 'bitscore': 2},
        top_score_percent=0.1,
        db_segment_size=10000, db_key='Taxonomy',
        ignore_ancestors=False, min_support=100, only_major=True,
        out_file='lca_analysis.txt',
        prefix_rank=True,
        show_path=False,
        list_reads=False)


def with_score():
    tre_file = 'resources/ncbi.tre'
    map_file = 'resources/ncbi.map'
    megan_map_file = 'resources/megan-map-Jan2021.db'
    blast_file = 'resources/Alice01-1mio-Jan-2021.txt'
    blast_map = {'qseqid': 0, 'sseqid': 1, 'bitscore': 2}
    top_score_percent = 0.1
    db_segment_size = 10000
    db_key = 'Taxonomy'
    ignore_ancestors = False
    min_support = 100
    only_major = False
    out_file = 'lca_analysis.txt'
    prefix_rank = True
    show_path = False
    list_reads = False

    # reads_ws, read_ids_ws = parse_blast_with_score(blast_file, blast_map)
    # reads_from_ws = filter_reads_by_top_score(reads_ws, 0.1)
    # reads, read_ids = parse_blast_filter(blast_file, 0.1, blast_map)

    # match = all([read_from_ws == read for read_from_ws, read in zip(reads_from_ws, reads)])
    # print('reads and reads with scores match: ' + str(match))

    # save_to_bin(reads, 'mapped_reads.bin')
    # save_to_bin(reads_ws, 'mapped_reads_ws.bin')

    # reads = map_accessions(reads, megan_map_file, db_segment_size, db_key)
    # reads_ws = map_accessions_with_scores(reads_ws, megan_map_file, db_segment_size, db_key)
    # reads_from_ws = filter_reads_by_top_score(reads_ws, top_score_percent)

    # save_to_bin(reads, 'mapped_reads.bin')
    # save_to_bin(reads_ws, 'mapped_reads_ws.bin')
    # save_to_bin(reads_from_ws, 'mapped_reads_from_ws.bin')

    # reads = load_from_bin('mapped_reads.bin')
    # reads_fws = load_from_bin('mapped_reads_from_ws.bin')


def project():
    tre_file = 'resources/ncbi.tre'
    map_file = 'resources/ncbi.map'
    megan_map_file = 'resources/megan-map-Jan2021.db'
    blast_file = 'resources/Alice01-1mio-Jan-2021.txt'
    blast_map = {'qseqid': 0, 'sseqid': 1, 'bitscore': 2}
    top_score_percent = 0.1
    db_segment_size = 10000
    db_key = 'Taxonomy'
    ignore_ancestors = False
    min_support = 200
    only_major = False
    out_file = 'lca_analysis.txt'
    prefix_rank = False
    show_path = False
    list_reads = True

    # mapped_reads = load_from_bin('mapped_reads.bin')
    # read_ids = load_from_bin('read_ids.bin')
    # tree = load_from_bin('tree.bin')
    # tree = parse_tree(tre_file, map_file)
    # id2address, address2id = compute_lca_addresses(tree)
    # map_lcas(tree, id2address, address2id, mapped_reads, read_ids, ignore_ancestors)
    tree = load_from_bin('tree_lca.bin')

    # project_reads_to_rank(tree, 'genus')
    apply_min_sup_filter(tree, min_support, only_major)

    # nodes = tree.nodes.values()
    # for node in nodes:
    #     if node.reads > 0:
    #         print(node.name, node.rank, node.reads)

    write_results(tree, out_file, prefix_rank, show_path, list_reads)


if __name__ == '__main__':
    project()

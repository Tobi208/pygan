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
        out_file='lca_analysis.txt')


def msf_only_major():
    run(tre_file='resources/ncbi.tre',
        map_file='resources/ncbi.map',
        megan_map_file='resources/megan-map-Jan2021.db',
        blast_file='resources/Alice01-1mio-Jan-2021.txt',
        blast_map={'qseqid': 0, 'sseqid': 1, 'bitscore': 2},
        top_score_percent=0.1,
        db_segment_size=10000, db_key='Taxonomy',
        ignore_ancestors=False, min_support=100, only_major=True,
        out_file='lca_analysis.txt')


def project():
    print('starting lca analysis')
    lca_start = time()
    tree = parse_tree('resources/ncbi.tre', 'resources/ncbi.map')
    id2address, address2id = compute_lca_addresses(tree)
    reads, read_ids =\
        parse_blast_filter('resources/Alice01-1mio-Jan-2021.txt', 0.1, {'qseqid': 0, 'sseqid': 1, 'bitscore': 2})
    mapped_reads = map_accessions(reads, 'resources/megan-map-Jan2021.db', 10000, 'Taxonomy')
    map_lcas(tree, id2address, address2id, mapped_reads, read_ids, False)
    project_reads_to_rank(tree, 'genus')
    write_results(tree, 'lca_analysis.txt')
    print('completed lca analysis in ' + timer(lca_start))


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


def pickler():

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

    reads_ws, read_ids_ws = parse_blast_with_score(blast_source, blast_map)
    reads, read_ids = parse_blast_filter(blast_source, 0.1, blast_map)

    save_to_bin(reads_ws, 'reads_ws.bin')
    reads_ws_bin = load_from_bin('reads_ws.bin')
    reads_from_ws_bin = filter_reads_by_top_score(reads_ws_bin, 0.1)

    match = all([read_from_ws_bin == read for read_from_ws_bin, read in zip(reads_from_ws_bin, reads)])
    print('reads and reads with scores from bin match: ' + str(match))


if __name__ == '__main__':
    with_score()
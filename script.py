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


def only_major():
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


if __name__ == '__main__':
    project()
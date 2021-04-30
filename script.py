from pygan.lca_analysis import *


def run_full_lca():
    run(tre_file='resources/ncbi.tre',
        map_file='resources/ncbi.map',
        megan_map_file='resources/megan-map-Jan2021.db',
        blast_file='resources/Alice01-1mio-Jan-2021.txt',
        blast_map={'qseqid': 0, 'sseqid': 1, 'bitscore': 2},
        top_score_percent=0.1,
        db_segment_size=10000, db_key='Taxonomy',
        ignore_ancestors=False, min_support=100,
        out_file='lca_analysis.txt')


if __name__ == '__main__':
    run_full_lca()
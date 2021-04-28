from pygan import lca_analysis

lca_analysis.execute(tre_file='resources/ncbi.tre',
                     map_file='resources/ncbi.map',
                     megan_map_file='resources/megan-map-Jan2021.db',
                     blast_file='resources/Alice01-1mio-Jan-2021.txt',
                     blast_format='tab', top_score_percent=0.1,
                     db_segment_size=10000, db_key='Taxonomy',
                     ignore_ancestors=False, min_support=100,
                     out_file='lca_analysis.txt')

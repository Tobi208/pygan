from pygan.tree.newick_parser import get_phylo_tree
from pygan.tree.map_parser import map_names
from pygan.blast.blast_parser import parse as blast_parse
from pygan.database.megan_map import fetch_all_taxonids


def lca_analysis():
    print('parsing tree')
    tree = get_phylo_tree('../resources/ncbi.tre')
    print('mapping names and ranks')
    map_names('../resources/ncbi.map', tree)
    print('parsing blast')
    all_accessions = blast_parse('../resources/Alice01-1mio-Jan-2021.txt', blast_format='tab', top_score_percent=0.1)
    print('fetching taxonids')
    all_taxonids = fetch_all_taxonids('../resources/megan-map-Jan2021.db', all_accessions)
    print('#reads: ' + str(len(all_taxonids)))


if __name__ == '__main__':
    lca_analysis()
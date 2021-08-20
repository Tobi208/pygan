# Python Metagenome Analzyer (PYGAN)

Requires Python 3.9 and a [MEGAN Map](https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html) to perform an LCA analysis on metagenomic sequences from alignment data. Refer to _LCA analysis of metagenomic sequences in Python_ or read the doc strings for further information.

## Run

Import the `pygan` module and call its `run` API to perfrom an analysis. 

```Python
pygan.run(tre_file='resources/ncbi.tre',
    map_file='resources/ncbi.map',
    megan_map_file='resources/megan-map-Jan2021.db',
    blast_file='resources/Alice01-1mio-Jan-2021.txt',
    blast_map={'qseqid': 0, 'sseqid': 1, 'bitscore': 2},
    top_score_percent=0.1,
    db_segment_size=10000, db_key='Taxonomy',
    ignore_ancestors=False, min_support=100, only_major=False,
    exclude=['genus'], project_mode='accession', project_rank='genus',
    cluster_degree=1, out_file='lca_analysis.txt',
    prefix_rank=True, show_path=False, list_reads=False)
```

### Description of the parameters

#### tre_file

Path to file containing taxonomy tree. NCBI or GTDB taxonomy recommended. 

#### map_file

Path to file containing mapping of taxonomy id to scientific name and rank. Must correspond to `tre_file`.

#### megan_map_file

Path to file containing [megan_map.db](https://software-ab.informatik.uni-tuebingen.de/download/megan6/welcome.html). It is highly recommended to store the database on a medium with fast reading (SSD).

#### blast_file

Path to file containing alignment data.

#### blast_map

Mapping of which column of the alignment data qseqid, sseqid and bitscore are in. E.g. `{'qseqid': 0, 'sseqid': 1, 'bitscore': 2}`. Required to parse the alignment data.

#### top_score_percent

Parameter used in the top score filter. Value must be between 0 and 1. An item within the percentage of the top score stays in the read. E.g. top score = 50, top_score_percent = 0.1: 47 remains, 43 is discarded. Prunes alignment data.

#### db_segment_size

Parameter to optimize performance of accession to taxon ID mapping. Different hardware may work better with different values. Recommended are values between 5,000 and 25,000.

#### db_key

Taxonomy to map accessions to. Use `'Taxonomy'` for NCBI and `'gtdb'` for GTDB.

#### ignore_ancestors

Whether to ignore ancestors in the LCA algorithm. When `True` no ancestor of another accesion can be the resulting lowest common ancestor.

#### min_support

Limit for the minimum support filter algorithm. All nodes with less reads than the minimum support limit forfeit their read to their parent nodes.

#### only_major

When `True` reads are forcibly pushed upwards to nodes of major ranks.

#### exclude

List of ranks to be exempted by the minimum support filter.

#### project_mode

Method of read projection to a specific rank. An attempt is made to contain all reads only in nodes of the target rank. Available methods are: ´'accession'´, ´'proportional'´, and ´'mixed'´ or leave empty `''`.

#### project_rank

Target rank for read projection.

#### cluster_degree

When choosing the `accession` method, specify cluster degree of at least 1. Potential hits are clustered. Reduces false-positives.

#### out_file

Path to output file. Output is generated as plain text.

#### prefix_rank

When `True` add an abbrevation of a node's rank to its name when printing.

#### show_path

When `True` print the entire path to the node instead of only its name

#### list_reads

When `True` print the list of a node's mapped read IDs. When `False` print the number of a node's mapped reads.


## Script

After familiarizing with the parameters and doc strings, script the analysis yourself or perform it in a REPL.

```Python
tree = parse_tree(tre_file, map_file)
id2address, address2id = compute_lca_addresses(tree)
reads, read_ids = parse_blast_filter(blast_file, top_score_percent, blast_map)
mapped_reads = map_accessions(reads, megan_map_file, db_segment_size, db_key)
map_lcas(tree, id2address, address2id, mapped_reads, read_ids, ignore_ancestors)
project_reads_to_rank(project_mode, tree, project_rank, mapped_reads, read_ids, cluster_degree)
apply_min_sup_filter(tree, min_support, exclude, only_major)
write_results(tree, out_file, prefix_rank, show_path, list_reads)
```

### Description of the methods

#### parse_tree

Takes the taxonomy data and produces an empty taxonomy tree.

#### compute_lca_addresses

Bidirectionally computes the address for each node in the phylogenetic tree. Later used in the LCA algorithm.

#### parse_blast_filter

Parse alignment data from a tabulated text file and apply the top score filter while parsing. To apply the top score filter specifically after parsing refer to `parse_blast_with_score`. Manually apply the top score filter with `filter_reads_by_top_score`. 

#### map_accessions

Maps accessions to taxons of the phylogenetic tree. Key must correspond to the specified taxonomy. Adjust `db_chunk_size` as necessary. Use `map_accessions_with_scores` to map reads that have not been filtered yet.

#### map_lcas

Populate the taxonomy with reads by applying the LCA algorithm.

#### project_reads

An attempt is made to contain all reads only in nodes of a specified rank. Available methods are: ´'accession'´, ´'proportional'´, and ´'mixed'´. 

#### apply_min_sup_filter

Applies the minimum support filter to the tree. Nodes that have less reads than the minimum support limit forfeit their reads to their parent nodes.

#### write_results

Generate plain text output of the taxonomy. Can prefix an abbrevation of the rank, show the entire path to the node and either list all read IDs or just show their number.


`save_to_bin` and `load_from_bin` allows for (de)serialization of data. May be useful to avoid multiple accession mappings or to store partial results of the analysis. Use `timer` to time your analysis duration.
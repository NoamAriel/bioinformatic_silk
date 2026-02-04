# Introduction

In the phylogenetic tree, there is an order called "Trichoptera" (caddisfly). 
Those insects living in the larva stage in fresh/creek waters, and in the adult stage they become moths.
My study focuses on them due to their special silk, which, in addition to the special properties of silk, including strength, extensibility, and toughness, is stable and adhesive in an aquatic environment.

Since there are many species of caddisfly, all producing silk with similar properties, it will be interesting to examine their amino acid composition and sequences to determine what features are essential for adhesive silk in an aquatic environment. 

For comparison, we will contrast these features with the silk of animals whose silk is not adhesive (negative control) and with those that produce adhesive silk (positive control).

Because available data on silk from aquatic animals is limited to a few species, each from a different order (Chironomus tepperi, Acentria ephemerella, and Argyroneta aquatica), the control will be based on terrestrial animals.

The analysis will involve examining the amino acid compositions of silks from species across different orders and searching for similar sequences and motifs among caddisfly species.
Those similarities are suspected to be related to their ability to produce adhesive silk in water. These findings will then be compared to sequences from other animals with known silk properties, such as silkworms, spiders, and bagworms.
Ultimately, this could help establish connections between sequences, amino acid composition, and the functional properties of the adhesive silk.


# Curernt code

## Order of operating

pls make sure you install the dependencies: ["requests","beautifulsoup4", "pandas", "numpy", "matplotlib"]

pls make sure you download the libraries folder

pls make sure you download the caddisfly, ants, moths and spiders folders.

pls make sure that each folder has the files: 

1.	ncbi_scrapper.py, 

2.	run_generate_species_index.py, 

3.	run_generate_taxonomy_graph.py 

4.	sxn_analysis_and_plotting.py.

This is also the operating order. 

pls notice that you update the path addresses correctly in each file.


## Code explanation 

The datasets of fibroin and spidroin proteins are downloaded by the ncbi_scrapper file from the NCBI database for each order. 

The files were characterized based on the protein type (such as heavy chain), and whether the file was reported as "partial" or not.

The taxonomy of each species is saved in MD and JSON files, named phylo_tree.

The run_generate_species_index file generates an index from all the download data

The run_generate_taxonomy_graph file generates a phylogenetic tree. 

The aa_composition_analysis file generates plots of the interesting amino acid composition, and tables that contain the data.

The sxn_analysis_and_plotting file generates plots of [SX]n motifs' percentages in the sequences, as well as X-residue compositions and tables containing the data.

The codes utilize a fixed color map for each amino acid, ensuring that the X-residue composition graph and future plots use the same color for each amino acid. That will help with an easier understanding of the different plots.

For convenient, the ncbi_scrapper, run_generate_species_index, run_generate_taxonomy_graph, aa_composition_analysis and sxn_analysis_and_plotting using libraries.

A standard deviation has been calculated for cases where proteins of a species have more than one data file, and it is presented as error bars in the plots.

There is a filtering system for the file run_generate_taxonomy_graph (phylogenetic tree generation), sxn_analysis_and_plotting, and aa_composition_analysis.

The filters are used to provide better control over presenting the data, resulting in more coherent, convenient, and understandable graphs. The specific filters are dependent on the graphs, since not all of them are relevant for all plots. 

## The filters
taxonomy_terms		 # filtering by the taxonomy name (trichoptera, for example). could be used by any taxonomy rank.  Default is [] / None

protein_types			# filtering by the type of protein (heavy chain, for example), default is [] / None. 

partial_full 			# filtering by sequence, if it is reported as partial or full. 

length_range 			# filtering by range sequence that interests, e.g., (100, 2450). Default is None.

length_threshold 		# filtering by threshold sequence that interests, e.g., 1500. default is None.

length_mode 			# determining if the threshold is “ge”: greater equal or “le”: less equal. The default is "ge". 

longest_factor = 2.0                  # optional default is 2.0, which means the shortest sequence can be at least half as long as the longest one.

longest_factor_scope  	# determining if the above length filters correspond to each species by itself, or all length sequences, or are relevant, even though they are from different species. In short, "species" (per organism) or "global" (all records). The default is "species".

## Special filters

In file” sxn_analysis_and_plotting”, in 273- 274 rows, there is an option to choose the range of n (min_n is the minimum motif length and max_n = 50 is the maximum motif length). There is not default. The user has to specify it.

In the file “run_generate_taxonomy_graph”, there is an option to choose the first and the last rank that the phylogenetic tree will present. Rank type (like superfamily) or ranks name (like Hydropsychoidea) are valid. The default is None.







MyTaxa
======

MyTaxa: an advanced taxonomy classifier for unknown genomic and metagenomic sequences

Author: Chengwei Luo, Luis Miguel Roderiguez-R, and Konstantinos Konstantinidis

Copyright (C): Konstantinidis Lab, Georgia Tech, 2013

Citation: Chengwei Luo, Luis Miguel Roderiguez-R, and Konstantinos T. Konstantinidis, MyTaxa: an advanced taxonomic classifier for genomic and metagenomic sequences. 2014, Nucleic Acids Research.

[What does it do?]
====================
MeTaxa represents a new algorithm that extends the Average Amino Acid Identity (AAI) concept (Konstantinidis and Tiedje, PNAS 2005) to identify the taxonomic affiliation of a query genome sequence or a sequence of a contig assembled from a metagenome, including short sequences (e.g., 100-1,000nt long), and to classify sequences representing novel taxa at three levels (whenever possible), i.e., species, genus and phylum. MeTaxa can assign a larger number of sequences and with higher accuracy compared to other tools available for the same purposes. This is largely attributed to the fact that MeTaxa considers all genes present in an unknown (query) sequence as classifiers and quantifies the classifying power of each gene using predetermined weights, which are derived from the analysis of orthologs of the gene from all available complete genomes. The weights are for i) how well the orthologs of the gene in question resolve the classification of the corresponding genomes at a given taxonomic level (species, genus, etc.) based on their degree of sequence conservation (for instance, the 16S rRNA gene resolves well at the genus and phylum levels but poorly at the species level); and ii) how frequently the ortholog gene phylogeny of the genomes compared deviates from the species phylogeny, the latter being approximated by the AAI tree, due primarily to horizontal gene transfer (HGT). MeTaxa also reports the statistical probability of the taxonomic assignment based on the Maximum Likelihood analysis. 

[Installation]
====================
cd to the MyTaxa directory in terminal, and do:

$ make

This will generate the executable binaries for 

If you haven't manually downloaded the pre-calculated database files and file them in /MyTaxa/db, you need to run:

$ python utils/download_db.py

and then decompress it by running:

$ tar xzvf db.latest.tar.gz

You should be all set after this.

[Usage]
====================
All you need is a gff file and a blast-like tabular file. In the gff file are the protein-coding gene annotations of the query sequences, and the tabular blast-like file has the results of searching against a reference database using those genes.

You then should run "utils/infile_convert.pl" to generate the input file for MyTaxa, and then you should run

$ MyTaxa [infile] [outfile] [thr] [num_hits]

thr is the threshold of scores (0-1) you define, and num_hits is the number of hits in the searching results to use (recommend 5)

The output is an XML style file with taxonomic information for each query sequence.

<strong>Please refer to the manual for detailed information on how to run it.</strong>

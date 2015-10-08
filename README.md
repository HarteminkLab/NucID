# NucID: Nucleosome Identification using DNase-seq
NucID is a Python package for identifying nucleosomes using single-end DNase-seq data:
  * It calculates a Bayes-factor-based nucleosome score for DNase-seq data in a 147 base pair genomic window. 
  * Genome-wide nucleosome scores are calculated using a sliding window approach. 
  * Nucleosomes are mapped using the genome-wide nucleosome score through a greedy algorithm.

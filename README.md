# NucID: Nucleosome Identification using DNase-seq
NucID is a Python package for identifying nucleosome positions from single-end DNase-seq data:
  * Using DNase-seq data, a Bayes-factor-based nucleosome score can be calculated for any 147 base pair genomic window. 
  * Genome-wide nucleosome scores are then determined by applying a simple sliding window approach. 
  * Finally, a greedy algorithm maps nucleosomes using these genome-wide nucleosome scores.

## Install
Clone this repository

```bash
git clone git@github.com:jianlingzhong/NucID.git
```

then

```bash
cd NucID && python setup.py install
```

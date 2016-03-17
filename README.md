# `NucID`: Nucleosome Identification using DNase-seq
`NucID` is a Python package for identifying nucleosome positions from single-end DNase-seq data:
  * Using DNase-seq data, a Bayes-factor-based nucleosome score can be calculated for any 147 base pair genomic window. 
  * Genome-wide nucleosome scores are then determined by applying a simple moving window approach. 
  * Finally, a greedy algorithm maps nucleosomes using these genome-wide nucleosome scores.

## Install
Clone this repository

```bash
git clone git@github.com:HarteminkLab/NucID.git
```

then

```bash
cd NucID && python setup.py install
```

## Usage
[example.ipynb](https://github.com/HarteminkLab/NucID/blob/master/example/example.ipynb) contains an example on how to calculate moving window nucleosome scores. Each file in the [package](https://github.com/HarteminkLab/NucID/tree/master/NucID) has an accompany `IPython` notebook that explains the code. 

## Pre-computed genome-wide nucleosome scores
We provide pre-computed genome-wide nucleosome scores for both *`S. cerevisiae`* genome (sacCer2) and *`human`* genome (hg18), available at:
  * [sacCer2](http://trackhub.genome.duke.edu/harteminklab/NucID/sacCer2/)
  * [hg18](http://trackhub.genome.duke.edu/harteminklab/NucID/hg18/)

Both were computed using the code here. 

## Reference
`NucID` was developed by [Jianling Zhong](https://github.com/jianlingzhong) in the Hartemink Lab for the following publication:

[Zhong, J., Luo, K., Winter, P. S., Crawford, G. E., Iversen E. S., & Hartemink, A. J. (2016). Mapping nucleosome positions using DNase-seq. Genome Res. gr.195602.115](http://genome.cshlp.org/content/early/2016/01/15/gr.195602.115.abstract)

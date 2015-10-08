
from distutils.core import setup

setup(
    name='NucID',
    version='0.1.0',
    description='Nucleosome Ientification using DNase-seq data: a Python package for mapping nucleosome translatinal position using DNase-seq data; an approach using Bayes factors',
    author='Jianling Zhong',
    url='http://jianlingzhong.github.io/NucID/',
    # package_data = {'': ['*.txt', '*.fasta', '*.gff', '*.csv']},
    packages=['NucID'],
    classifiers=[
        'Intended Audience :: Academic',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Cython',
        'Topic :: Scientific/Engineering :: Statistics',
    ],

)

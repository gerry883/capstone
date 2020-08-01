Uniprot proteomes for sequence autocorrelation descriptors

Input files:
1. [proteome.fasta.gz] the organism proteome fasta.gz file of the uniprot database (https://www.uniprot.org/proteomes/?query=&sort=score)

Output file:
1. [descriptors.csv] the csv file containing PyBiomed Autocorrelation descriptors

Required software:
1. Autocorrelation.py, PyProtein.py
2. conda environment with python 3, numpy, pandas, scikitlearn
3. fastaparser (conda install -c kronopt fastaparser)

Usage:
1. python generate.py proteome.fasta.gz descriptors.csv



Kaggle dataset to csv

Intructions:
1. create the env
conda create --name capstone python=2.7 matplotlib scikit-learn
conda install -c conda-forge rdkit (or conda install -c rdkit rdkit on linux)
conda install -c openbabel openbabel
2. download the https://github.com/gadsbyfly/PyBioMed/raw/master/PyBioMed/download/PyBioMed-1.0.zip
3. extract or uncompress the PyBioMed-1.0.zip file
4. cd PyBioMed-1.0
python setup.py install
5. download the dataset here https://www.kaggle.com/shahir/protein-data-set/download
    
Input files:
1. [11797_16251_bundle_archive.zip] the zip containing the two dataset

Output file:
1. [pr_classification.csv] the csv file containing PyBiomed Autocorrelation descriptors

Usage:
1. python kaggle_to_csv.py 11797_16251_bundle_archive.zip pr_classification.csv

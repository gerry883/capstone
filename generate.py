import os
import sys
import numpy as np
from pipeline import *

fasta_proteome_path = sys.argv[1]
outpath = sys.argv[2]

data = pipeline(fasta_proteome_path)
# outfile = open(outpath, 'w')
# data.to_csv(outfile)
# outfile.close()
data.to_csv(outpath)
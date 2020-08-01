import pandas as pd
import numpy as np
# from PyBioMed.PyGetMol import GetProtein
from PyProtein import PyProtein
import gzip
import fastaparser
limit = None



def fasta_parser(path):
    fasta_df = pd.DataFrame()
    print("Done")
    with gzip.open(path,'rt') as fasta_file:
        reader = fastaparser.Reader(fasta_file)
        print("File loaded")
        fasta_df['ID'] = [seq.id for seq in reader]
        print("Column 1 of 3 added")
        fasta_df['Description'] = [seq.description for seq in reader]
        print("Column 2 of 3 added")
        fasta_df['Sequence'] = [seq.sequence_as_string() for seq in reader]
        print("Column 3 of 3 added")
        
    return fasta_df

def clean_from_fasta_final(fasta_df):
    
    aa = ["A", "R", "N", "D", "C", "E",
           "Q", "G", "H", "I", "L", "K",
           "M", "F", "P", "S", "T", "W",
           "Y", "V"]
    
    # process the df
    # fasta_df.drop('Unnamed: 0',axis=1,inplace=True)
    mask = fasta_df.Sequence.apply(lambda x: set(x) <= set(aa))
    data = fasta_df[mask]
    print('Fasta additional amino-acids filtered out')
    
    return data


def to_seq(data):
    
    sequences = data.Sequence.tolist()
    indexes = data.ID
    
    return sequences, indexes


def create_array(sequences, indexes, limit=limit):
    
    # to generate the dataframe columns
    protein = PyProtein(sequences[0])
    #indexes = indices.to_numpy(copy=True)
    columns = []    
    columns.extend(protein.GetMoreauBrotoAuto().keys())
    columns.extend(protein.GetMoranAuto().keys())
    columns.extend(protein.GetGearyAuto().keys())
    
    # to create the dataframe rows
    rows = indexes[:limit]
    
    # array created with rows-keys shape
    array = np.zeros([len(rows),len(columns)])
    
    # descriptor generation
    for row in range(len(rows)):
        protein_class = PyProtein(sequences[row])
        dicti = {}
        dicti.update(protein_class.GetMoreauBrotoAuto())
        dicti.update(protein_class.GetMoranAuto())
        dicti.update(protein_class.GetGearyAuto())
        
        # array value assignation
        for key in dicti.keys():
            array[row, columns.index(key)] = dicti[key]
            
        # report
        if row%100 == 0:
            print('Done '+str(row) + ' of '+ str(len(rows)))
        elif row == len(rows):
            print('Done '+ str(row))
            
    return array, rows, columns
    
    
def descr_dataframe_autocorrelation(array, rows, columns):
    
    ## to dataframe
    df = pd.DataFrame(array, columns=columns, index=rows)
    
    return df

def pipeline(path):
    
    fasta_df = fasta_parser(path)
    print("Done")
    data = clean_from_fasta_final(fasta_df)
    print('Done')
    sequences, indexes = to_seq(data)
    print('Done')
    array, rows, columns = create_array(sequences, indexes)
    print('Done')
    df = descr_dataframe_autocorrelation(array, rows, columns)
    print('Done')
    
    return df
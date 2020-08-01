import pandas as pd
import numpy as np
from PyBioMed.PyGetMol import GetProtein
from PyBioMed import Pyprotein
from zipfile import ZipFile
import os
import sys
limit = None

    
def descr_dataframe_autocorrelation(array, rows, columns):
    
    ## to dataframe
    df = pd.DataFrame(array, columns=columns, index=rows)
    
    return df

def create_array(sequences, indexes, limit=limit):
    
    # to generate the dataframe columns
    protein = Pyprotein.PyProtein(sequences[0])
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
        protein_class = Pyprotein.PyProtein(sequences[row])
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


def clean2(zippath, filenames):
    
    # amino-acid letters
    aa = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    
    # files to dataframe
    with ZipFile(zippath, 'r').open(filenames[0]) as file1:
        df_seq = pd.read_csv(file1)        
    with ZipFile(zippath, 'r').open(filenames[1]) as file2:
        df_dup = pd.read_csv(file2)
    
    # merge the dataframes
    df_merge = df_dup.merge(df_seq,how='inner',on='structureId')
    df_merge.rename({'macromoleculeType_x':'macromoleculeType',
                                            'residueCount_y':'residueCount'},axis=1,inplace=True)
    
    # clean the columns
    df_merge.drop(['macromoleculeType_y','residueCount_x'],axis=1,inplace=True)
    df_merge_protein = df_merge[df_merge.macromoleculeType=='Protein']
    df_merge_proteins = df_merge_protein[df_merge_protein.chainId=='A']
    df_merge_proteins2 = df_merge_proteins.drop([
                         'experimentalTechnique',
                         'chainId',
                         'residueCount',
                         'macromoleculeType',
                         'crystallizationMethod',
                         'crystallizationTempK',
                         'densityMatthews',
                         'densityPercentSol',
                         'pdbxDetails',
                         'phValue',
                         'publicationYear',
                         'resolution'], axis=1)
    
    # drop na
    df_merge_proteins3 = df_merge_proteins2.dropna()
    
    # to identify the top 10 classes
    tab = df_merge_proteins2.classification.value_counts(normalize=True)
    filtered_classes = list(tab.keys()[:10])
    df_merge_proteins3 = df_merge_proteins2[df_merge_proteins2.classification.isin(filtered_classes)]
    
    # to delete unreadable amino-acids letter
    mask = df_merge_proteins3.sequence.apply(lambda x: set(x) <= set(aa))
    df_masked = df_merge_proteins3[mask]
    
    # cleansed sequences to list
    sequences = df_masked.sequence.tolist()
    
    # indexes
    indexes = df_masked.structureId
    return sequences, indexes


def pipeline(zip_file, file_list):
    
    sequences, indexes = clean2(zip_file, file_list)
    array, rows, columns = create_array(sequences, indexes)
    df = descr_dataframe_autocorrelation(array, rows, columns)
    return df
    
zip_file = sys.argv[1]
outfile = sys.argv[2]
file_list = ['pdb_data_seq.csv','pdb_data_no_dups.csv']

df = pipeline(zip_file, file_list)
df.to_csv(outfile)



import numpy as np
#import numba as nb
import anndata as an
import time
#from combat.pycombat import pycombat
import sklearn as sk
import pandas as pd
import scipy.io       as sc_io
import scipy.sparse   as sc_sp
import scipy.optimize as sc_opt

SEARCH_DATA_PATH = '/Users/andrewleduc/Desktop/nPOP_Paper/dat/raw_data/pSCoPE/ev_updated.txt'
META_DATA_PATH = '/Users/andrewleduc/Desktop/nPOP_Paper/dat/raw_data/pSCoPE/annotation.csv'
CELL_TYPES = ['m','u']
NEG_CTRL = ['neg']

def preproc_DIAN(data_path):
    search = pd.read_csv(data_path, sep='\t')#,engine= 'pyarrow'

def preproc_MQ(data_path1,data_path2):
    #List of columns for Reporter ion quant
    RI_col = [str(x) for x in range(2, 19)]
    RI_col = ['Reporter intensity ' + s2 for s2 in RI_col]

    #Additional important features to read in
    var = ['Raw file','Modified sequence','Leading razor protein',
           'Charge','PIF','PEP','dart_qval','Reverse','Potential contaminant']

    #Read in raw and meta data
    search_dat = pd.read_csv(data_path1, sep='\t',low_memory=False,usecols=var+RI_col)#,engine= 'pyarrow'
    search_dat.columns = search_dat.columns.str.replace(' ', '_')
    meta = pd.read_csv(data_path2, sep=',',low_memory=True)#,engine= 'pyarrow'

    #Data filtering
    search_dat = search_dat[search_dat.Raw_file.isin(meta['Set'].unique())]
    search_dat = search_dat[search_dat.PIF > .8 ]
    search_dat = search_dat[search_dat.dart_qval < .01]
    search_dat = search_dat[search_dat.Potential_contaminant != '+']
    search_dat = search_dat[search_dat.Reverse != '+']

    #reference normalization
    search_dat.iloc[:, 6:23] = search_dat.iloc[:, 6:23].div(search_dat.Reporter_intensity_2, axis=0)

    #Condense data to matrix form
    var = [element.replace(' ', '_') for element in var]
    RI_col = [element.replace(' ', '_') for element in RI_col]
    search_dat = pd.melt(search_dat, id_vars= var , value_vars= RI_col)
    search_dat['Sample'] = search_dat['Raw_file']+search_dat['variable']
    new_cols = var + ['Sample'] + ['value']
    search_dat = search_dat[new_cols]
    search_dat = search_dat.pivot_table(index=['Leading_razor_protein','Modified_sequence'], columns='Sample', values = 'value')

    #Map meta data to MQ results
    meta = pd.melt(meta, id_vars= ['Set'] , value_vars= RI_col)
    meta['Sample'] = meta['Set'] + meta['variable']
    meta = meta[meta.Sample.isin(search_dat.columns)]
    id_vect = [str(x) for x in range(1, (len(meta.index)+1))]
    id_vect = ['i' + x for x in id_vect]
    meta['id'] = id_vect
    meta = meta[meta.value.isin(CELL_TYPES+NEG_CTRL)]

    search_dat = search_dat[meta['Sample']]
    search_dat.columns = meta['id']

    return search_dat,meta

def filter_cell_CV(df,meta,cutoff):

    df_sc = df[meta['id']]
    df_sc = df_sc.to_numpy()
    df_sc = df_sc / np.median(df_sc, axis=0)
    df_sc = df_sc / np.mean(df_sc, axis=1)

    CV_mat = np.zeros(len(df['Leading_razor_protein'].unique()),shape(df_sc)[1])
    prot_list_uni = df['Leading_razor_protein'].unique()
    prot_list = np.array(df['Leading_razor_protein'])

    count = 0
    for prot in prot_list_uni:

        arr_mask = np.where(prot_list == prot)
        arr_index = np.arange(0, len(prot_list))[arr_mask]

        df_temp = df_sc[arr_index][:]
        CV_mat[count][:] = np.std(df_temp, axis=0)/np.mean(df_temp, axis=0)

        count += 1

    sc_keep = np.median(CV_mat, axis=0) < cutoff

    return sc_keep

def normalize_mat(df,sc_ids,Protein_quant):
    df = df[sc_ids]
    df_sc = df.to_numpy()

    # Filter for peptides/ single cells with correct amount of missing data
    df_sc = df_sc[np.sum(df_sc,axis=1)/len(df_sc.index) > .95][:]
    df_sc = df_sc[:][np.sum(df_sc, axis=0)/len(df_sc.columns) > .95]

    # Keep meta data updated
    df = df[]

    # Normalize peptide data
    df = df/np.median(df,axis=0)
    df = df / np.mean(df, axis=1)

    #Collapse to protein level
    if Protein_quant == 1:
        sdf

    # Normalize protein data
    df = df / np.median(df, axis=0)
    df = df / np.mean(df, axis=1)

    #Impute
    imputer = sk.impute.KNNImputer(n_neighbors=3)
    imputer.fit_transform(df)

    #Batch correct
    #data_corrected = pycombat(data, batch)

def first_pass_PCA():
    print('hi')
    #print Rsq due to different batch effects first 5 PCs
    #PCA no imp/with imp

#change to anndata

if __name__ == '__main__':

    #matrix= np.arange(0,9).reshape((3, 3))

    #print(np.std(matrix, axis=0) / np.mean(matrix, axis=0))

    '''
    start_time = time.time()
    df1,meta = preproc_MQ(SEARCH_DATA_PATH,META_DATA_PATH)
    sc_keep = filter_cell_CV(df1, meta, .4)
    
    
    print("--- %s seconds ---" % (time.time() - start_time))
    '''




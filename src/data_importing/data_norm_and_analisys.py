import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from inmoose.pycombat import pycombat_norm
import os
from sklearn.preprocessing import RobustScaler
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
from helpers import get_first_indexs, apply_KNN_impute,box_plot


def run_preprocessing(
    plot_nan = False,
    plot_boxPlots = True,
    no_change = False):

    path = PROCESSED_DATA_FOLDER

    out_path = FILTERING_FIGURES

    # Create directory if it doesn't exist
    os.makedirs(out_path, exist_ok=True)

    try:
        filtered_df = pd.read_csv(path+'filter.csv', index_col=0)
        print('succesfully loaded data')
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        big_df = pd.read_csv(COMBINED_DATA_OUTPUT_FILE,index_col=0)
        # filter the data on 20%
        nan_percentage_genes = big_df.isna().mean(axis=1) * 100

        # Filter rows where NaN percentage <= 20%
        filtered_genes = nan_percentage_genes[nan_percentage_genes <= 20].index

        # Keep only the filtered rows
        filtered_df = big_df.loc[filtered_genes]


        # Calculate the percentage of NaN values in each column
        nan_percentage_samples = filtered_df.isna().mean() * 100

        # Filter columns where NaN percentage <= 20%
        filtered_columns = nan_percentage_samples[nan_percentage_samples <= 20].index

        # Keep only the filtered columns
        filtered_df = filtered_df[filtered_columns]

        filtered_df.to_csv(path+'filter.csv')

    print('data loaded')
    big_df = filtered_df

    #! Remove columns from the datafram
    #TODO: static
    # big_df.loc[:, ~(big_df > 1000000).any()]
    big_df = big_df.filter(regex=r'^(?!.*GSM463683).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463684).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463685).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463686).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463687).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463688).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463689).*$')
    big_df = big_df.filter(regex=r'^(?!.*GSM463690).*$')

    matrix = big_df.to_numpy()
    matrix_nan = big_df.isna().to_numpy()

    chromosomes = ['1','2','3','4','5']
    indices:list[int] = get_first_indexs(big_df.index,chromosomes)

    if plot_nan:
        row_nan_count = big_df.isna().sum(axis=1)
        col_nan_count = big_df.isna().sum(axis=0)
        
        plt.imshow(matrix_nan, cmap='hot', interpolation='nearest')
        plt.savefig(out_path+'matrix.svg')
        plt.close()
        # Create directory if it doesn't exist
        os.makedirs(out_path+'row', exist_ok=True)
        os.makedirs(out_path+'col', exist_ok=True)
        for i,c in enumerate(indices):
            min_var = indices[i]
            try:
                max_var = indices[i+1]
            except:
                max_var = len(matrix)
            plt.bar(range(len(row_nan_count.index[min_var:max_var])),row_nan_count.values[min_var:max_var])
            plt.ylim(0,1850)
            plt.xlabel('Genes')
            plt.ylabel('Missing number of data points')
            plt.savefig(out_path+'row/0.row_dis'+chromosomes[i]+'.svg')
            plt.close()

        plt.bar(range(len(col_nan_count.index)),col_nan_count.values)
        plt.xlabel('Experiments')
        plt.ylabel('Missing number of data points')
        plt.savefig(out_path+'col/0.col_dis.svg')
        plt.close()


    # With simple replacement
    print('plotted Nans', plot_nan)
    np.nan_to_num(matrix,copy=False)
    print('nans filled with 0')

    for column in big_df:
        # print(big_df[column])
        if (big_df[column]>300).any():
            big_df[column] = np.log2(big_df[column]+1)

    print('applied log2')
    def get_study(sample: str):
        return sample.split('_')[0]
    study_map = list(map(get_study,big_df.columns))

    # KNN Impute
    try: 
        print('reading file from: ' + path+'imputed.csv')
        df_impute:pd.DataFrame = pd.read_csv(path+'/imputed.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        print('Could not read file, running KNN impute')
        df_impute:pd.DataFrame = apply_KNN_impute(big_df,5)
        print('KNN impute ran, saving file')
        df_impute.to_csv(path+'imputed.csv')
        print('file saved at: ' + path+'imputed.csv')
        # get the UMAP


    # NORMALIZE
    #? Apply batch correction
    try:
        study_corrected_df:pd.DataFrame = pd.read_csv(path+'study_corrected.csv',index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        study_map = list(map(get_study,df_impute.columns))
        # raise ValueError("Running in the cluster")
        d = dict([(y,x+1) for x,y in enumerate(sorted(set(study_map)))])
        batches = []
        for el in study_map:
            batches.append(d[el])
        study_corrected_df:pd.DataFrame = pycombat_norm(df_impute, batches)
        study_corrected_df.to_csv(path+'study_corrected.csv')
    try:
        standardized_df = pd.read_csv(path+'/standardized.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        standardized_df = ((df_impute.T - df_impute.T.mean()) / df_impute.T.std()).T
        standardized_df.to_csv(path+'/standardized.csv')
    try:
        robust_df = pd.read_csv(path+'/robust.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError
        scaler = RobustScaler()
        robust_df = pd.DataFrame(scaler.fit_transform(df_impute), columns=df_impute.columns)
        robust_df.to_csv(path+'/robust.csv')
    
    try:
        double_norm = pd.read_csv(path+'/2_way_norm.csv', index_col=0)
    except FileNotFoundError:
        if no_change:
            raise FileNotFoundError            
        mat = study_corrected_df.to_numpy()
        q75, q25 = np.percentile(mat, [75 ,25],axis=1,keepdims=True)
        iqr = q75 - q25
        norm = (mat - np.median(mat,axis=1,keepdims=True))/iqr
        double_norm = pd.DataFrame(norm, columns=study_corrected_df.columns,index=study_corrected_df.index)
        double_norm.to_csv(path+'/2_way_norm.csv')

    
    if plot_boxPlots:
        box_plot(df_impute,120, out_path+'/uncorrected_box_plots/')
        box_plot(study_corrected_df,120, out_path+'/study_corrected_box_plots/')
        box_plot(double_norm,120, out_path+'/2_way_norm_plots/')
        box_plot(robust_df,120, out_path+'/robust/')
        box_plot(standardized_df,120,out_path+'/standardized/')
    print('Done')

if __name__ == '__main__':
    run_preprocessing()
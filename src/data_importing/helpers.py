
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import umap
import os
import math
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.impute import KNNImputer
import matplotlib.cm as cm
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn
from matplotlib.colors import LogNorm
import sys

def get_geo_list(path:str):
    read =  pd.read_csv(path)
    read = read.loc[read['depository_source'] == 'GEO']
    read = read.loc[read['species'] == 'Arabidopsis thaliana']
    return list(read['depository_accession'])
def mapping(x):
    if type(x) is str:
        return x.upper()
    else:
        return x
    
def predicate(gene:str, chromosome:str)-> bool:
    return str('AT'+chromosome+'G') in gene

def get_first_indexs(df_index,chromo:list[str]):
    array = []
    for i in chromo:
        gene:str = next(filter(lambda x : predicate(x,str(i)), df_index))
        array.append(df_index.get_loc(gene))
    return array


def get_Umap(matrix: np.array, study_map: list = None, name: str = '',
             save_loc: str = '', title: str = 'UMAP projection of the dataset'):
    """
    Generate and save UMAP projection plot, creating directories if needed.
    
    Args:
        matrix: Input data matrix
        study_map: List of labels for coloring points (optional)
        name: Additional identifier for output filename
        save_loc: Directory to save the plot
        title: Plot title
    """
    # Create directory if it doesn't exist
    os.makedirs(save_loc, exist_ok=True)
    
    # Perform UMAP transformation
    reducer = umap.UMAP()
    scaled_data = StandardScaler().fit_transform(matrix)
    embedding = reducer.fit_transform(scaled_data)
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    if study_map is None:
        plt.scatter(
            embedding[:, 0],
            embedding[:, 1]
        )
    else:
        colors = cm.rainbow(np.linspace(0, 1, max(study_map)+1))
        for num, emb in enumerate(embedding):
            plt.scatter(
                emb[0],
                emb[1],
                color=colors[study_map[num]]
            )
    
    plt.gca().set_aspect('equal', 'datalim')
    plt.title(title, fontsize=24)
    
    # Construct save path and save figure
    output_path = os.path.join(save_loc, f'umap{name}.svg')
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()
    
    # print(f"UMAP plot saved to: {output_path}")


def normalize(arr, t_min, t_max):
    norm_arr = []
    diff = t_max - t_min
    diff_arr = max(arr) - min(arr)    
    for i in arr:
        temp = (((i - min(arr))*diff)/diff_arr) + t_min
        norm_arr.append(temp)
    return norm_arr

def normalize_2d(matrix):
    norm = np.linalg.norm(matrix)
    matrix = matrix/norm  # normalized matrix
    return matrix

def plot_sim_matrix(matrix: np.array, indices: list = None, chromosomes: list = None, 
                   name: str = '', save_loc: str = '', title: str = ''):
    """
    Plot similarity matrix and save to specified location, creating directories if needed.
    
    Args:
        matrix: Input data matrix
        indices: List of indices to split the matrix
        chromosomes: List of chromosome names for labeling
        name: Additional name identifier for output file
        save_loc: Base directory to save outputs
        title: Plot title
    """
    # Determine folder structure
    folder = 'Genes/'
    if indices is None:
        indices = [0]
        folder = 'Samples/'

    if chromosomes is None:
        chromosomes = ['']

    # Create directories if they don't exist
    output_dir = os.path.join(save_loc, 'sim_matrix', folder)
    os.makedirs(output_dir, exist_ok=True)

    for i, _ in enumerate(indices):
        # print('Plotting sim matrix', i)
        min_idx = indices[i]
        try:
            max_idx = indices[i+1]
        except IndexError:
            max_idx = len(matrix)
        
        # print('Computing similarity')
        # Compute pairwise cosine similarity
        similarity_matrix = cosine_similarity(matrix[min_idx:max_idx])
        
        # print('Creating plot')
        plt.imshow(similarity_matrix, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.title(title)
        
        # Construct output path
        output_path = os.path.join(output_dir, f'sim_{chromosomes[i]}_matrix{name}.svg')
        plt.savefig(output_path)
        plt.close()
        # print(f'Finished plot saved to {output_path}')

    plt.close()
    # print('Done with all similarity plots')

def plot_heat_map(df:pd.DataFrame,save_loc:str, name: str,cluster:bool=True):
    # Create directories if they don't exist
    output_dir = os.path.join(save_loc, 'heat_map')
    os.makedirs(output_dir, exist_ok=True)
    o = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    seaborn.clustermap(df,row_cluster=cluster, col_cluster=False,method='complete', norm=LogNorm())
    plt.savefig(output_dir+'/'+name+'.png')
    plt.close()
    sys.setrecursionlimit(o)

def apply_KNN_impute(df:pd.DataFrame,n_neighbors: int):
    imputer = KNNImputer(n_neighbors=n_neighbors)

    # Fit and transform the dataset
    df_imputed = pd.DataFrame(imputer.fit_transform(df), columns=df.columns, index=df.index)
    return df_imputed


def hierarchical_clustering(data_matrix:np.array):
    linkage_data = linkage(data_matrix, method='ward', metric='euclidean')
    return linkage_data

def hierarchical_clustering_plot(data_matrix:np.array, path:str, name:str):
    linkage_data = linkage(data_matrix, method='ward', metric='euclidean', optimal_ordering=True)
    dendrogram(linkage_data, no_labels= True)
    plt.savefig(path+'cluster_ordered_'+name+'.svg')
    plt.close()

def box_plot(df: pd.DataFrame, cols_per_plot: int, out_path: str):
    """
    Generates and saves boxplots from a DataFrame, with visual separators between studies.
    
    Args:
        df (pd.DataFrame): The input data.
        cols_per_plot (int): The number of columns (samples) to include in each plot.
        out_path (str): The directory to save the output plots.
    """
    num_cols = len(df.columns)
    num_plots = math.ceil(num_cols / cols_per_plot)
    
    # Ensure the output directory exists
    os.makedirs(out_path, exist_ok=True)

    for plot_num in range(num_plots):
        start_idx = plot_num * cols_per_plot
        end_idx = min((plot_num + 1) * cols_per_plot, num_cols)
        
        current_cols = df.iloc[:, start_idx:end_idx]
        
        # Dynamically adjust figure width based on the number of columns
        # This prevents the plots from looking too compressed.
        fig_width = max(20, len(current_cols.columns) * 0.5) 
        plt.figure(figsize=(fig_width, 10))
        
        # Set a fixed y-axis limit for consistent comparison across plots
        plt.ylim(-18, 18)
        
        # Create the boxplot
        plt.boxplot(current_cols, labels=current_cols.columns)
        
        # --- NEW: Add vertical lines between studies ---
        # Extract study IDs from column names (e.g., 'GSE12345' from 'GSE12345_sample1')
        study_ids = [name.split('_')[0] for name in current_cols.columns]
        
        # Iterate through the columns to find where the study ID changes
        for i in range(1, len(study_ids)):
            if study_ids[i] != study_ids[i-1]:
                # Add a vertical line. Positions are 1-based, so the line goes
                # at i + 0.5 to be between box i and box i+1.
                plt.axvline(x=i + 0.5, color='black', linestyle='--', linewidth=1)
        # --- END NEW ---

        # Rotate x-axis labels for better readability
        plt.xticks(rotation=45, ha='right')
        
        # Add title and adjust layout to prevent labels from being cut off
        plt.title(f'Boxplot Group {plot_num + 1} (Columns {start_idx + 1}-{end_idx})')
        plt.tight_layout()
        
        # Save the figure and close it to free up memory
        plt.savefig(os.path.join(out_path, f'boxplot_group_{plot_num + 1}.png'))
        plt.close()
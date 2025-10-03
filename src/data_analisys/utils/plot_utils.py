import matplotlib.pyplot as plt
import torch
import networkx as nx
import os
import sys
import seaborn
from matplotlib.colors import LogNorm
import pandas as pd
import umap
import numpy as np
from sklearn.manifold import TSNE
from scipy.cluster import hierarchy
from typing import Optional

def plot_losses(losses,output_dir,exp_name,iteration):
    plt.close()
    xs = [x for x in range(len(losses))]
    plt.plot(xs, losses)
    plt.savefig(f'{output_dir}{exp_name}_loss_{iteration}.svg')
    plt.close()

def plot_values(values,output_dir,exp_name):
    xs = [x for x in range(len(values))]
    plt.plot(xs, values)
    plt.savefig(f'{output_dir}/{exp_name}.svg')
    plt.close()

def plot_values_bar(values,output_dir,exp_name,title:str='',y_label='',x_label=''):
    xs = [x for x in range(len(values))]
    plt.bar(xs, values)
    plt.title(title)
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    plt.savefig(f'{output_dir}/{exp_name}.svg')
    plt.close()


def plot_heat_map(df:pd.DataFrame,save_loc:str, name: str,cluster:bool=True,typ:str='png',title='',log_norm:bool = True, col: list = None,col_cluster:bool=False):
    # Create directories if they don't exist
    output_dir = os.path.join(save_loc, 'heat_map')
    os.makedirs(output_dir, exist_ok=True)
    o = sys.getrecursionlimit()
    sys.setrecursionlimit(10000)
    test = seaborn.clustermap(df,row_cluster=cluster, col_cluster=col_cluster,method='complete', norm=LogNorm() if log_norm else None, col_colors=col)
    # plt.title(title, fontsize=24)
    plt.savefig(f'{output_dir}/{name}.{typ}')
    plt.close()
    sys.setrecursionlimit(o)
    return test

def plot_predictions(final_predictions,final_targets,N,output_dir,exp_name,x_name='',y_name=''):
    xs_ = [x for x in range(len(final_predictions[0]))]
    for plo in range(len(final_predictions)):
        torch.save((final_predictions[plo]/N), output_dir+f'saves/{plo}_pred.pt')
        torch.save((final_targets[plo]/N), output_dir+f'saves/{plo}_targ.pt')
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.plot(xs_, (final_predictions[plo]/N),label='prediction')
        plt.plot(xs_, (final_targets[plo]/N), label = 'target')
        plt.legend(loc="upper left")
        plt.savefig(output_dir+exp_name+'_'+str(plo)+'.svg')
        plt.close()
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.plot(xs_, (final_predictions[plo]/N),label='prediction')
        plt.legend(loc="upper left")
        plt.savefig(output_dir+'pred_'+exp_name+'_'+str(plo)+'.svg')
        plt.close()

def plot_weights(plot_mat,output_dir,exp_name,x_name='',y_name=''):
    for i,w in enumerate(plot_mat):
        m = [x for x in range(len(w))]
        plt.xlabel(x_name)
        plt.ylabel(y_name)
        plt.plot(m, w.cpu(),label='weights')
        plt.legend(loc="upper left")
        plt.savefig(output_dir+exp_name+'_weights_'+str(i)+'.svg')
        plt.close()
    # plt.savefig(output_dir+exp_name+'_weights_overlap.svg')
    # plt.close()

def plot_matrix_graph(plot_mat_bool,output_dir,exp_name,remove_isolated:bool = False,x_name='',y_name=''):
    G = nx.from_numpy_array(plot_mat_bool.numpy())
    if remove_isolated:
        isolated_nodes = list(nx.isolates(G))
        G.remove_nodes_from(isolated_nodes)
    nx.draw(G,node_size=30, alpha = 0.5)
    plt.axis('equal')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.savefig(output_dir+exp_name+'_graph.svg')
    plt.close()

def plot_matrix(mat,location,name,x_name='',y_name='',title=''):
    plt.imshow(mat, cmap='hot', interpolation='nearest')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.title(title)
    plt.colorbar()
    plt.savefig(f'{location}{name}.svg')
    plt.close()


def get_Umap(matrix: np.array, colour_map: Optional[list] = None, name: str = '',
             save_loc: str = './outputs/misc', title: str = 'UMAP projection of the dataset', small_dots=False):
    """
    Generate and save UMAP projection plot, creating directories if needed.
    
    Args:
        matrix: Input data matrix
        colour_map: List of labels for coloring points (optional)
        name: Additional identifier for output filename
        save_loc: Directory to save the plot
        title: Plot title
        small_dots: Boolean flag for small dot size
    """
    # Create directory if it doesn't exist
    os.makedirs(save_loc, exist_ok=True)
    
    # Perform UMAP transformation
    reducer = umap.UMAP(n_neighbors=10)
    embedding = reducer.fit_transform(matrix)
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Define marker and color options
    markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', '8', 'P', 'd', '>', '<', '1', '2', '3', '4']
    colors = plt.cm.rainbow(np.linspace(0, 1, len(set(colour_map)) if colour_map is not None else ['b']))
    
    if colour_map is None:
        plt.scatter(
            embedding[:, 0],
            embedding[:, 1],
            s=0.5 if small_dots else 10,
            c='b',
            marker='o'
        )
    else:
        unique_groups = sorted(set(colour_map))
        for group in unique_groups:
            # Get indices of samples in this group
            group_indices = [i for i, x in enumerate(colour_map) if x == group]
            
            # Assign marker and color
            marker_idx = group % len(markers)
            color_idx = group % len(colors)
            
            plt.scatter(
                embedding[group_indices, 0],
                embedding[group_indices, 1],
                c=[colors[color_idx]] * len(group_indices),
                marker=markers[marker_idx],
                s=0.5 if small_dots else 10,
                label=f'Group {group}'
            )
        
        # Add legend if there are groups
        if len(unique_groups) > 1:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.gca().set_aspect('equal', 'datalim')
    plt.title(title, fontsize=24)
    
    # Construct save path and save figure
    output_path = os.path.join(save_loc, f'umap_{name}.svg')
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()

def get_Umap_3(embedding, 
             colour_map: Optional[list] = None, 
             marker_map: Optional[list] = None,
             name: str = '',
             save_loc: str = './outputs/misc', 
             title: str = 'UMAP projection of the dataset', 
             small_dots: bool = False,
             legend: bool = True):
    """
    Generate and save UMAP projection plot with independent color and marker assignments.
    
    Args:
        embedding: Umap embedings
        colour_map: List of color labels for points
        marker_map: List of marker labels for points
        name: Additional identifier for output filename
        save_loc: Directory to save the plot
        title: Plot title
        small_dots: Boolean flag for small dot size
        legend: Whether to show legend
    """
    # Create directory if it doesn't exist
    os.makedirs(save_loc, exist_ok=True)
    
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Define marker and color options
    available_markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', '8', 'P', 'd', '>', '<', '1', '2', '3', '4']
    available_colors = plt.cm.rainbow(np.linspace(0, 1, len(set(colour_map)) if colour_map is not None else ['b']))
    
    if colour_map is None and marker_map is None:
        # Simple scatter plot if no mappings provided
        plt.scatter(
            embedding[:, 0],
            embedding[:, 1],
            s=0.5 if small_dots else 15,
            c='b',
            marker='o'
        )
    else:
        # Create combined groups based on unique color-marker pairs
        if colour_map is None:
            colour_map = [0] * len(embedding)
        if marker_map is None:
            marker_map = [0] * len(embedding)
            
        unique_combinations = sorted(set(zip(colour_map, marker_map)))
        
        # Create mappings for colors and markers
        color_dict = {color: available_colors[i % len(available_colors)] 
                     for i, color in enumerate(sorted(set(colour_map)))}
        marker_dict = {marker: available_markers[i % len(available_markers)] 
                       for i, marker in enumerate(sorted(set(marker_map)))}
        
        # Plot each unique combination
        for color_val, marker_val in unique_combinations:
            mask = (np.array(colour_map) == color_val) & (np.array(marker_map) == marker_val)
            plt.scatter(
                embedding[mask, 0],
                embedding[mask, 1],
                c=[color_dict[color_val]] * sum(mask),
                marker=marker_dict[marker_val],
                s=0.5 if small_dots else 15,
                label=f'Color {color_val}, Marker {marker_val}' if legend else None
            )
    plt.gca().set_aspect('equal', 'datalim')
    plt.title(title, fontsize=24)
    #TODO: plot the colours of the dots to the legend
    if True and (colour_map is not None or marker_map is not None):
        # Create custom legend handles with both color and marker information
        handles = []
        if colour_map is not None and marker_map is not None:
            # For both color and marker mappings
            uni = np.unique(np.array(colour_map),return_counts = True)
            for color_val, marker_val in unique_combinations:
                handles.append(plt.Line2D([0], [0], 
                                        marker=marker_dict[marker_val], 
                                        color='w',
                                        markerfacecolor=color_dict[color_val],
                                        markersize=10,
                                        label=f'{color_val} ({marker_val}) #{uni[1][list(uni[0]).index(color_val)]}'))
        elif colour_map is not None:
            # For color mapping only
            for color_val in sorted(set(colour_map)):
                handles.append(plt.Line2D([0], [0], 
                                        marker='o', 
                                        color='w',
                                        markerfacecolor=color_dict[color_val],
                                        markersize=10,
                                        label=f'{color_val}'))
        elif marker_map is not None:
            # For marker mapping only
            for marker_val in sorted(set(marker_map)):
                handles.append(plt.Line2D([0], [0], 
                                        marker=marker_dict[marker_val], 
                                        color='w',
                                        markerfacecolor='gray',
                                        markersize=10,
                                        label=f'{marker_val}'))
        
        plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    
    # Construct save path and save figure
    output_path = os.path.join(save_loc, f'umap_{name}.svg')
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()
    return embedding

def get_Umap_2(matrix: np.array, 
             colour_map: Optional[list] = None, 
             marker_map: Optional[list] = None,
             name: str = '',
             save_loc: str = './outputs/misc', 
             title: str = 'UMAP projection of the dataset', 
             small_dots: bool = False,
             legend: bool = True):
    """
    Generate and save UMAP projection plot with independent color and marker assignments.
    
    Args:
        matrix: Input data matrix
        colour_map: List of color labels for points
        marker_map: List of marker labels for points
        name: Additional identifier for output filename
        save_loc: Directory to save the plot
        title: Plot title
        small_dots: Boolean flag for small dot size
        legend: Whether to show legend
    """
    # Create directory if it doesn't exist
    os.makedirs(save_loc, exist_ok=True)
    
    # Perform UMAP transformation
    reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5)
    embedding = reducer.fit_transform(matrix)
    
    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Define marker and color options
    available_markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', '8', 'P', 'd', '>', '<', '1', '2', '3', '4']
    available_colors = plt.cm.rainbow(np.linspace(0, 1, len(set(colour_map)) if colour_map is not None else ['b']))
    
    if colour_map is None and marker_map is None:
        # Simple scatter plot if no mappings provided
        plt.scatter(
            embedding[:, 0],
            embedding[:, 1],
            s=0.5 if small_dots else 15,
            c='b',
            marker='o'
        )
    else:
        # Create combined groups based on unique color-marker pairs
        if colour_map is None:
            colour_map = [0] * len(embedding)
        if marker_map is None:
            marker_map = [0] * len(embedding)
            
        unique_combinations = sorted(set(zip(colour_map, marker_map)))
        
        # Create mappings for colors and markers
        color_dict = {color: available_colors[i % len(available_colors)] 
                     for i, color in enumerate(sorted(set(colour_map)))}
        marker_dict = {marker: available_markers[i % len(available_markers)] 
                       for i, marker in enumerate(sorted(set(marker_map)))}
        
        # Plot each unique combination
        for color_val, marker_val in unique_combinations:
            mask = (np.array(colour_map) == color_val) & (np.array(marker_map) == marker_val)
            plt.scatter(
                embedding[mask, 0],
                embedding[mask, 1],
                c=[color_dict[color_val]] * sum(mask),
                marker=marker_dict[marker_val],
                s=0.5 if small_dots else 15,
                label=f'Color {color_val}, Marker {marker_val}' if legend else None
            )
    plt.gca().set_aspect('equal', 'datalim')
    plt.title(title, fontsize=24)
    #TODO: plot the colours of the dots to the legend
    if True and (colour_map is not None or marker_map is not None):
        # Create custom legend handles with both color and marker information
        handles = []
        if colour_map is not None and marker_map is not None:
            # For both color and marker mappings
            uni = np.unique(np.array(colour_map),return_counts = True)
            for color_val, marker_val in unique_combinations:
                handles.append(plt.Line2D([0], [0], 
                                        marker=marker_dict[marker_val], 
                                        color='w',
                                        markerfacecolor=color_dict[color_val],
                                        markersize=10,
                                        label=f'{color_val} ({marker_val}) #{uni[1][list(uni[0]).index(color_val)]}'))
        elif colour_map is not None:
            # For color mapping only
            for color_val in sorted(set(colour_map)):
                handles.append(plt.Line2D([0], [0], 
                                        marker='o', 
                                        color='w',
                                        markerfacecolor=color_dict[color_val],
                                        markersize=10,
                                        label=f'{color_val}'))
        elif marker_map is not None:
            # For marker mapping only
            for marker_val in sorted(set(marker_map)):
                handles.append(plt.Line2D([0], [0], 
                                        marker=marker_dict[marker_val], 
                                        color='w',
                                        markerfacecolor='gray',
                                        markersize=10,
                                        label=f'{marker_val}'))
        
        plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
    
    # Construct save path and save figure
    output_path = os.path.join(save_loc, f'umap_{name}.svg')
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()
    return embedding



def plot_tsne(
    df: pd.DataFrame,
    colors: list,
    markers:list,
    title: str = "t-SNE Projection",
    name: str = '',
    point_size: int = 30,
    alpha: float = 0.7,
    legend: bool = False,
    save_path: Optional[str] = None,
):
    """
    Plot t-SNE visualization with custom colors for each data point.

    Args:
        df: Input DataFrame (samples x features).
        colors: List of colors (hex/rgb/names) for each data point.
        title: Plot title.
        figsize: Figure dimensions.
        point_size: Size of scatter points.
        alpha: Transparency of points.
        legend: Whether to show a legend (if colors map to categories).
        save_path: If provided, saves the plot to this path.
    """
    os.makedirs(save_path, exist_ok=True)
    # Compute t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    tsne_results = tsne.fit_transform(df.values)

    # Create plot
    plt.figure(figsize=(10, 8))
    
    # Define marker and color options
    available_markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', '8', 'P', 'd', '>', '<', '1', '2', '3', '4']
    available_colors = plt.cm.rainbow(np.linspace(0, 1, len(set(colors)) if colors is not None else ['b']))
    
    if colors is None and markers is None:
        # Simple scatter plot if no mappings provided
        plt.scatter(
            tsne_results[:, 0],
            tsne_results[:, 1],
            c='b',
            marker='o'
        )
    else:
        # Create combined groups based on unique color-marker pairs
        if colors is None:
            colors = [0] * len(tsne_results)
        if markers is None:
            colors = [0] * len(tsne_results)
            
        unique_combinations = sorted(set(zip(colors, markers)))
        
        # Create mappings for colors and markers
        color_dict = {color: available_colors[i % len(available_colors)] 
                     for i, color in enumerate(sorted(set(colors)))}
        marker_dict = {marker: available_markers[i % len(available_markers)] 
                       for i, marker in enumerate(sorted(set(markers)))}
        
        # Plot each unique combination
        for color_val, marker_val in unique_combinations:
            mask = (np.array(colors) == color_val) & (np.array(markers) == marker_val)
            plt.scatter(
                tsne_results[mask, 0],
                tsne_results[mask, 1],
                c=[color_dict[color_val]] * sum(mask),
                marker=marker_dict[marker_val],
                label=f'Color {color_val}, Marker {marker_val}' if legend else None
            )
    
    plt.gca().set_aspect('equal', 'datalim')
    plt.title(title, fontsize=24)

    
    # Construct save path and save figure
    plt.savefig(f'{save_path}/{name}.svg', format='svg', bbox_inches='tight')
    plt.close()


def plot_dendogram(embedding,linkage_method,number_of_clusters,figure_out_path,name=''):
        os.makedirs(f'{figure_out_path}/dendogram/', exist_ok=True)
        # 1. Compute linkage matrix
        Z = hierarchy.linkage(embedding, method=linkage_method)

        # 2. Plot dendrogram
        plt.figure(figsize=(12, 8))
        dendro = hierarchy.dendrogram(
            Z,
            truncate_mode='lastp',  # show only the last p merged clusters
            p=number_of_clusters,   # show only these many clusters
            show_leaf_counts=True,  # show number of samples in each cluster
            leaf_rotation=90.,      # rotate labels for better readability
            leaf_font_size=12.,     # font size for labels
            show_contracted=True,   # show contracted branches
            color_threshold=Z[-number_of_clusters+1, 2]  # color threshold for clusters
        )

        plt.title(f'{name} Dendrogram ({linkage_method} linkage)')
        plt.xlabel('Sample index or (cluster size)')
        plt.ylabel('Distance')
        plt.grid(False)
        plt.tight_layout()
        plt.axhline(y=Z[-number_of_clusters+1, 2], color='r', linestyle='--')
        plt.savefig(f'{figure_out_path}/dendogram/dendogram_{name}.svg')
        plt.close()

def evaluate_cluters(clusters,out_path):
    _,counts = np.unique(clusters,return_counts=True)
    plt.plot(counts)
    plt.savefig(f'{out_path}/test.svg')
    return
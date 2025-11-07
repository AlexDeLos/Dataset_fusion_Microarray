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
import re

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


def plot_heat_map(df:pd.DataFrame,save_loc:str, name: str,cluster:bool=True,typ:str='png',title='',log_norm:bool = True, col: Optional[pd.DataFrame] = None,col_cluster:bool=False):
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

def plot_summary_scores(scores_dict: dict, title: str, file_name: str, output_dir: str):
    """
    Generates an improved bar plot with grouped section labels for clustering scores.
    Labels like 'robust val TREATMENT' are split into a main group label 'robust'
    and a bar label 'val TREATMENT'.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'bar_rand_ind.svg').
        output_dir (str): The directory to save the plot in.
    """
    if not scores_dict:
        print(f"Skipping plot '{title}' because the score dictionary is empty.")
        return

    # --- 1. Prepare Data and Labels ---
    full_labels = list(scores_dict.keys())
    values = list(scores_dict.values())

    # Split labels into the prefix (group) and the rest (bar label)
    # e.g., 'robust val TREATMENT' -> group='robust', short_label='val TREATMENT'
    processing_types = [label.split()[0] for label in full_labels]
    short_labels = [' '.join(label.split()[1:]) for label in full_labels]

    # --- 2. Define Colors ---
    color_map = {
        'TREATMENT': '#1f77b4', 'TISSUE': '#2ca02c', 'study': '#d62728',
    }
    # Color bars based on the last word of the original full label
    bar_colors = [color_map.get(label.split()[-1], '#7f7f7f') for label in full_labels]

    # --- 3. Create Plot ---
    # Use subplots for better control over the layout
    fig, ax = plt.subplots(figsize=(16, 10))
    ax.bar(range(len(values)), values, align='center', color=bar_colors)

    # --- 4. Add Dotted Separators and Section Labels ---
    # Find where the processing type changes to draw lines and place labels
    change_indices = [i for i in range(1, len(processing_types)) if processing_types[i] != processing_types[i-1]]
    
    # Define the start and end of each section
    section_starts = [0] + change_indices
    section_ends = change_indices + [len(full_labels)]

    # Increase the bottom margin to make space for the new section labels
    fig.subplots_adjust(bottom=0.25)

    for start, end in zip(section_starts, section_ends):
        # Draw a vertical line at the start of each new section (except the first)
        if start > 0:
            ax.axvline(x=start - 0.5, color='gray', linestyle='--', linewidth=1.5)

        # Calculate the midpoint of the section to place the group label
        mid_point = start + (end - 1 - start) / 2
        label_text = processing_types[start]
        
        # mapping to show paper friendly labels
        name_map ={
            'robust': 'Robust Normalization',
            'standardized': 'Standardized',
            'robust+': 'Robust Norm two waysalization',
            'standardized+': 'Standardized two ways',
            '2_way_norm': 'Study corrected then genes Robust Normalizedalization',
            'study_corrected': 'Study corrected',
            'imputed': 'No correction'
        }
        label_text = name_map[label_text]
        # Add the group label text below the x-axis
        ax.text(mid_point, -0.35, label_text, ha='center', va='top', fontsize=14, weight='bold',
                transform=ax.get_xaxis_transform(),rotation = 45)

    # --- 5. Final Plot Adjustments ---
    ax.set_xticks(range(len(full_labels)))
    ax.set_xticklabels(short_labels, rotation=90, fontsize=10)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title(title, fontsize=16)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Ensure the output directory exists and save the figure
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, file_name))
    plt.close()

def plot_iteration_scores(scores_dict: dict, y_label: str, title: str, file_name: str, output_path: str):
    """
    Generates and saves a colored bar plot for scores from a single processing iteration.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        y_label (str): The label for the Y-axis.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'silhouette_scores.svg').
        output_path (str): The directory path to save the plot in.
    """
    if not scores_dict:
        print(f"Skipping plot '{title}' because the score dictionary is empty.")
        return

    labels = list(scores_dict.keys())
    values = list(scores_dict.values())

    # Define a color map for different label types
    color_map = {
        'TREATMENT': '#1f77b4',  # Muted blue
        'TISSUE': '#2ca02c',     # Green
        'study': '#d62728',      # Red
    }
    # Assign a color to each bar based on its label, with a default gray
    bar_colors = [color_map.get(label, '#7f7f7f') for label in labels]

    # --- Create the plot ---
    plt.figure(figsize=(10, 7))
    plt.bar(labels, values, color=bar_colors)

    # --- Final plot adjustments ---
    plt.ylabel(y_label, fontsize=12)
    plt.title(title, fontsize=16)
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, file_name))
    plt.close()


def plot_iteration_scores_modified(scores_dict: dict, y_label: str, title: str, file_name: str, output_path: str):
    """
    Generates and saves a colored bar plot for scores, restructured to compare
    Treatment, Tissue, and Study scores for each matrix type.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        y_label (str): The label for the Y-axis.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'rand_index_scores.svg').
        output_path (str): The directory path to save the plot in.
    """
    if not scores_dict:
        print(f"Skipping plot '{title}' because the score dictionary is empty.")
        return

    # 1. Data Restructuring and Parsing
    parsed_data = {}
    for label, score in scores_dict.items():
        # Regex to capture the three parts: Matrix_Val, Classification Target
        # The pattern looks for: (matrix_name val) (Classification Target)
        # Note: study scores usually don't have 'val' but this regex covers it
        match = re.search(r'(.+?(?:val|study)) (\w+)$', label)
        if match:
            matrix_name_val = match.group(1).strip() # e.g., 'robust val' or 'robust study'
            target = match.group(2).strip().lower()  # e.g., 'treatment' or 'tissue'
        else:
            # Handle cases that don't fit the expected pattern gracefully (e.g., 'imputed val study')
            print(f"Warning: Could not fully parse label '{label}'. Skipping...")
            continue
        
        # We need a unique key for each group of 3 bars (e.g., 'robust val')
        group_key = matrix_name_val.replace(' val', '').replace(' study', '')
        
        if group_key not in parsed_data:
            parsed_data[group_key] = {'treatment': 0, 'tissue': 0, 'study': 0, 'order': len(parsed_data)}

        # Map the original target to the new, simplified key
        simplified_target = target
        if simplified_target not in ['treatment', 'tissue', 'study']:
             # Fallback for unexpected targets, though likely not needed here
             simplified_target = 'study'

        parsed_data[group_key][simplified_target] = score

    # 2. Prepare Data for Plotting
    # Sort groups based on their appearance in the original dictionary
    sorted_groups = sorted(parsed_data.items(), key=lambda item: item[1]['order'])
    
    # Create lists for plotting
    matrix_names = []  # The labels to appear once under the group of bars
    scores = []        # The actual values for all bars
    colors = []        # The color for each bar
    
    # Define an explicit color map for the simplified targets
    color_map_simplified = {
        'treatment': '#1f77b4',  # Muted blue
        'tissue': '#2ca02c',     # Green
        'study': '#d62728',      # Gray (for comparison)
    }

    # Order the scores for plotting: Treatment, Tissue, Study for each matrix
    for matrix_name, data in sorted_groups:
        matrix_names.append(matrix_name) # Used for tick labels

        # Append scores and colors in the desired order
        scores.extend([data['treatment'], data['tissue'], data['study']])
        colors.extend([color_map_simplified['treatment'], 
                       color_map_simplified['tissue'], 
                       color_map_simplified['study']])
        
    # The X-axis labels for individual bars will just be 'treatment', 'tissue', 'study'
    x_bar_labels = (['treatment', 'tissue', 'study'] * len(matrix_names))
    
    # 3. Create the Plot
    
    plt.figure(figsize=(12, 7)) # Increased size for better readability
    
    # Generate positions for the bars
    bar_width = 0.25
    x_positions = []
    
    # Calculate positions for grouped bars
    for i, _ in enumerate(matrix_names):
        center = i * (3 * bar_width + 0.2) # Group separation
        x_positions.append(center - bar_width)     # Treatment position
        x_positions.append(center)                 # Tissue position
        x_positions.append(center + bar_width)     # Study position

    plt.bar(x_positions, scores, color=colors, width=bar_width, align='center')

    # --- Final plot adjustments ---
    plt.ylabel(y_label, fontsize=12)
    plt.title(title, fontsize=16)
    
    # Set the main tick labels (the matrix names) to be centered under the group
    main_tick_positions = [i * (3 * bar_width + 0.2) for i, _ in enumerate(matrix_names)]
    plt.xticks(main_tick_positions, matrix_names, rotation=45, ha="right", fontsize=10)

    # Add a second, minor x-axis for 'treatment', 'tissue', 'study' labels
    ax = plt.gca()
    # Create the minor ticks (one for each bar)
    ax.set_xticks(x_positions, minor=True)
    # Set the minor tick labels (the bar types)
    ax.set_xticklabels(['TRT', 'TIS', 'STD'] * len(matrix_names), 
                       minor=True, rotation=90, ha='center', fontsize=8) 
    
    # Minor tick labels need to be adjusted to appear *under* the main labels
    # This requires a bit of manual tweaking often, but for this structure, we'll
    # rely on the minor ticks property for a basic look.
    
    # Add grid lines
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    
    # Add vertical lines to delineate the matrix groups clearly
    for i in range(1, len(matrix_names)):
         plt.axvline(x=main_tick_positions[i] - (1.5 * bar_width), color='gray', linestyle='--', alpha=0.4)
    
    plt.tight_layout()
    
    # Ensure the output directory exists
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        
    plt.savefig(os.path.join(output_path, file_name))
    plt.close()


import matplotlib.pyplot as plt
import os
import re

def plot_summary_scores_modified(scores_dict: dict, title: str, file_name: str, output_dir: str):
    """
    Generates an improved bar plot with grouped section labels for clustering scores,
    with bars re-ordered to group TREATMENT, TISSUE, and study scores consecutively
    for each matrix type.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'bar_rand_ind_grouped.svg').
        output_dir (str): The directory to save the plot in.
    """
    if not scores_dict:
        print(f"Skipping plot '{title}' because the score dictionary is empty.")
        return

    # --- 1. Data Restructuring and Parsing ---
    parsed_data = {}
    original_order_keys = []
    
    # Regex to capture the three parts: Matrix_Name, Validation_Status, Classification Target
    # Example: '2_way_norm_og val TREATMENT'
    # Group 1: (.+?(?:val|study)) -> '2_way_norm_og val' or 'robust study'
    # Group 2: (\w+)$ -> 'TREATMENT', 'TISSUE', or 'study'
    
    # Pre-parse and structure the data
    for label, score in scores_dict.items():
        match = re.search(r'(.+?(?:val)) (\w+)$', label)
        
        # Special case: handling 'study' scores which often don't have 'val'
        if not match:
             match = re.search(r'(.+?) (study)$', label)
        
        if not match:
            print(f"Warning: Could not fully parse label '{label}'. Skipping...")
            continue
            
        # Example: 'robust val' or 'robust'
        matrix_val = match.group(1).strip()
        # Example: 'treatment' or 'tissue' or 'study'
        target = match.group(2).strip().lower()

        # Determine the group key (Matrix Name) and ensure it's recorded only once
        group_key = matrix_val.replace(' val', '').replace(' val study', '').replace(' study', '')
        
        if group_key not in parsed_data:
            # Initialize with default values and record insertion order
            parsed_data[group_key] = {'treatment': 0, 'tissue': 0, 'study': 0, 'order': len(original_order_keys)}
            original_order_keys.append(group_key)

        # Map the original target to the simplified key (treatment, tissue, study)
        simplified_target = target
        if simplified_target == 'treatment':
            parsed_data[group_key]['treatment'] = score
        elif simplified_target == 'tissue':
            parsed_data[group_key]['tissue'] = score
        elif simplified_target == 'study':
            parsed_data[group_key]['study'] = score
        else:
             # This should ideally not happen if targets are only TREATMENT/TISSUE/study
             print(f"Warning: Unhandled target '{target}' for group '{group_key}'.")


    # 2. Re-order and Flatten Data for Plotting
    
    full_labels = []   # The simplified labels for the x-axis (TREATMENT, TISSUE, study)
    values = []        # The actual scores
    bar_colors = []    # The color for each bar
    
    # Define an explicit color map for the simplified targets
    color_map = {
        'treatment': '#1f77b4',  # Blue
        'tissue': '#2ca02c',     # Green
        'study': '#d62728',      # Gray
    }
    
    # List to track the group type for the section labels
    processing_types = []
    
    # Sort groups based on their original appearance
    for group_key in original_order_keys:
        data = parsed_data[group_key]
        
        # Plot the three targets consecutively: Treatment, Tissue, Study
        values.extend([data['treatment'], data['tissue'], data['study']])
        full_labels.extend(['TREATMENT', 'TISSUE', 'study'])
        bar_colors.extend([color_map['treatment'], color_map['tissue'], color_map['study']])
        
        # Repeat the group key three times for section label tracking
        processing_types.extend([group_key, group_key, group_key])
        
    # --- 3. Create Plot ---
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.bar(range(len(values)), values, align='center', color=bar_colors)

    # --- 4. Add Dotted Separators and Section Labels ---
    
    # Find where the processing type changes (which is every 3 bars in the new structure)
    change_indices = [i for i in range(1, len(processing_types)) if processing_types[i] != processing_types[i-1]]
    
    section_starts = [0] + change_indices
    section_ends = change_indices + [len(full_labels)]

    # Increase the bottom margin to make space for the new section labels
    fig.subplots_adjust(bottom=0.25)

    # mapping to show paper friendly labels for the group names
    name_map ={
        'robust': 'Robust Norm',
        'standardized': 'Standardized',
        'robust+': 'Robust Norm two ways',
        'standardized+': 'Standardized two ways',
        '2_way_norm_og': 'Study corrected then genes Robust Normalized OG',
        '2_way_norm': 'Study corrected then genes Robust Normalized',
        'study_corrected': 'Study corrected',
        'imputed': 'No correction'
    }

    for start, end in zip(section_starts, section_ends):
        # Draw a vertical line at the start of each new section (except the first)
        if start > 0:
            ax.axvline(x=start - 0.5, color='gray', linestyle='--', linewidth=1.5)

        # Calculate the midpoint of the section to place the group label
        mid_point = start + (end - 1 - start) / 2
        
        # The key is the same for all three bars in the section
        raw_label_text = processing_types[start] 
        label_text = name_map.get(raw_label_text, raw_label_text)
        
        # Add the group label text below the x-axis
        ax.text(mid_point, -0.35, label_text, ha='center', va='top', fontsize=12, weight='bold',
                transform=ax.get_xaxis_transform(), rotation=45)

    # --- 5. Final Plot Adjustments ---
    ax.set_xticks(range(len(full_labels)))
    # Set the simplified bar labels
    ax.set_xticklabels(full_labels, rotation=90, fontsize=9)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title(title, fontsize=16)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Ensure the output directory exists and save the figure
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, file_name),bbox_inches='tight')
    plt.close()


def plot_summary_scores_relative(scores_dict: dict, title: str, file_name: str, output_dir: str):
    """
    Generates a bar plot showing scores for TREATMENT and TISSUE relative to the
    Study score (X - Study Score). Study scores are not plotted.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'bar_rand_ind_relative.svg').
        output_dir (str): The directory to save the plot in.
    """
    if not scores_dict:
        print(f"Skipping plot '{title}' because the score dictionary is empty.")
        return

    # --- 1. Data Restructuring, Parsing, and Calculation ---
    parsed_data = {}
    original_order_keys = []
    
    # First pass: Gather all scores including Study to calculate the difference
    for label, score in scores_dict.items():
        # Capture Matrix_Name and Classification Target
        match = re.search(r'(.+?(?:val|study)) (\w+)$', label)
        
        if not match:
            print(f"Warning: Could not fully parse label '{label}'. Skipping...")
            continue
            
        # Example: 'robust val' or 'imputed val study'
        matrix_val = match.group(1).strip()
        target = match.group(2).strip().lower()

        # Determine the group key (Matrix Name)
        # Handle cases like 'imputed val study' where 'study' is part of the name
        group_key = matrix_val.replace(' val', '').replace(' val study', '').replace(' study', '')
        
        if group_key not in parsed_data:
            # Initialize with default values and record insertion order
            parsed_data[group_key] = {'treatment': 0, 'tissue': 0, 'study': 0, 'order': len(original_order_keys)}
            original_order_keys.append(group_key)

        # Store the raw score
        parsed_data[group_key][target] = score


    # Second pass: Calculate relative scores (Score - Study Score)
    for group_key, data in parsed_data.items():
        study_score = data['study']
        
        # Calculate Treatment and Tissue scores relative to Study
        data['treatment'] = data['treatment'] - study_score
        data['tissue'] = data['tissue'] - study_score

    # --- 2. Flatten Data for Plotting (Excluding Study) ---
    
    full_labels = []   # The simplified labels for the x-axis (TREATMENT, TISSUE)
    values = []        # The relative scores (X - Study Score)
    bar_colors = []    # The color for each bar
    
    # Define an explicit color map 
    color_map = {
        'treatment': '#1f77b4',  # Blue
        'tissue': '#2ca02c',     # Green
    }
    
    # List to track the group type for the section labels
    processing_types = []
    
    # Sort groups based on their original appearance
    for group_key in original_order_keys:
        data = parsed_data[group_key]
        
        # Plot the two targets consecutively: Treatment, Tissue
        values.extend([data['treatment'], data['tissue']])
        full_labels.extend(['TREATMENT', 'TISSUE'])
        bar_colors.extend([color_map['treatment'], color_map['tissue']])
        
        # Repeat the group key twice for section label tracking
        processing_types.extend([group_key, group_key])
        
    # --- 3. Create Plot ---
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.bar(range(len(values)), values, align='center', color=bar_colors)

    # Add a horizontal line at y=0 for reference, as scores are now relative
    ax.axhline(0, color='black', linewidth=0.8, linestyle='-')

    # --- 4. Add Dotted Separators and Section Labels ---
    
    # Find where the processing type changes (which is every 2 bars in the new structure)
    change_indices = [i for i in range(1, len(processing_types)) if processing_types[i] != processing_types[i-1]]
    
    section_starts = [0] + change_indices
    section_ends = change_indices + [len(full_labels)]

    # Increase the bottom margin to make space for the new section labels
    fig.subplots_adjust(bottom=0.25)

    # mapping to show paper friendly labels for the group names
    name_map ={
        'robust': 'Robust Norm',
        'standardized': 'Standardized',
        'robust+': 'Robust Norm two ways',
        'standardized+': 'Standardized two ways',
        '2_way_norm_og': 'Study corrected then genes Robust Normalized OG',
        '2_way_norm': 'Study corrected then genes Robust Normalized',
        'study_corrected': 'Study corrected',
        'imputed': 'No correction'
    }

    for start, end in zip(section_starts, section_ends):
        # Draw a vertical line at the start of each new section (except the first)
        if start > 0:
            ax.axvline(x=start - 0.5, color='gray', linestyle='--', linewidth=1.5)

        # Calculate the midpoint of the section to place the group label
        mid_point = start + (end - 1 - start) / 2
        
        raw_label_text = processing_types[start] 
        label_text = name_map.get(raw_label_text, raw_label_text)
        
        # Add the group label text below the x-axis
        ax.text(mid_point, -0.35, label_text, ha='center', va='top', fontsize=12, weight='bold',
                transform=ax.get_xaxis_transform(), rotation=45)

    # --- 5. Final Plot Adjustments ---
    ax.set_xticks(range(len(full_labels)))
    ax.set_xticklabels(full_labels, rotation=90, fontsize=9)
    # Update Y-label to reflect the relative score
    ax.set_ylabel(f'{ax.get_ylabel().replace("Score", "Relative Score (X - Study Score)")}', fontsize=12)
    ax.set_title(title, fontsize=16)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Ensure the output directory exists and save the figure
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, file_name),bbox_inches='tight')
    plt.close()
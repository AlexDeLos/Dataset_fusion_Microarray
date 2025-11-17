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

def plot_projection(embedding,
    colors: list,
    markers: list,
    title: Optional[str] = "t-SNE Projection",
    name: str = '',
    legend: bool = True,
    save_path: str = ''):

    os.makedirs(save_path, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 10))
    # Define marker and color options
    available_markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'X', 'h', '8', 'P', 'd', '>', '<', '1', '2', '3', '4']
    
    # Handle None cases for labels
    if colors is None:
        colors = ['default'] * len(tsne_results)
    if markers is None:
        markers = ['default'] * len(tsne_results)
        
    unique_colors = sorted(set(colors))
    unique_markers = sorted(set(markers))
    
    available_colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_colors)))

    # Create mappings for colors and markers
    color_dict = {color: available_colors[i % len(available_colors)] 
                  for i, color in enumerate(unique_colors)}
    marker_dict = {marker: available_markers[i % len(available_markers)] 
                   for i, marker in enumerate(unique_markers)}
    
    # Create combined groups based on unique color-marker pairs
    unique_combinations = sorted(set(zip(colors, markers)))
    
    # Plot each unique combination
    for color_val, marker_val in unique_combinations:
        # Handle the default case
        if color_val == 'default' and marker_val == 'default':
            mask = np.ones(len(embedding), dtype=bool)
            plot_color = 'b'
            plot_marker = 'o'
            plot_label = None
        else:
            mask = (np.array(colors) == color_val) & (np.array(markers) == marker_val)
            plot_color = color_dict[color_val]
            plot_marker = marker_dict[marker_val]
            # Your original code only labels by color, we'll keep that logic
            plot_label = f'{color_val}' if legend else None

        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            c=[plot_color] * sum(mask),
            marker=plot_marker,
            label=plot_label
        )
    
    ax.set_aspect('auto', 'datalim','C')
    if title:
        ax.set_title(title, fontsize=24)

    if legend:
        # 1. Get unique handles and labels
        # This de-duplicates the labels (e.g., 'Tissue A' appearing 5 times)
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles)) # de-duplicate
        
        # 2. Automatically determine number of columns
        # Aim for max 20 rows per column (adjust as needed)
        num_items = len(unique)
        max_rows_per_col = 40  
        num_cols = (num_items + max_rows_per_col - 1) // max_rows_per_col
        num_cols = max(1, num_cols) # Ensure at least 1 column

        # 3. Place legend outside the plot
        ax.legend(
            unique.values(), 
            unique.keys(),
            loc='upper left',
            bbox_to_anchor=(1.02, 1.0), # (102% width, 100% height)
            borderaxespad=0.0,
            ncol=num_cols,
            fontsize='medium' # 'small' or 'medium' is good for papers
        )
    # plt.axis('square')
    # Construct save path and save figure
    if save_path:
        # We MUST use bbox_inches='tight' so the legend isn't cut off
        fig.savefig(f'{save_path}/{name}.svg', format='svg', bbox_inches='tight')
    
    plt.close(fig) # Close the figure object

def plot_tsne(
    df: pd.DataFrame,
    colors: list,
    markers: list,
    save_path: str,
    title: Optional[str] = "t-SNE Projection",
    name: str = '',
    legend: bool = True):
    """
    Plot t-SNE visualization with an improved legend placed outside the plot.

    Args:
        df: Input DataFrame (samples x features).
        colors: List of labels to be used for color mapping.
        markers: List of labels to be used for marker mapping.
        title: Plot title.
        name: Filename for saving.
        legend: Whether to show a legend.
        save_path: If provided, saves the plot to this path.
    """        
    # Compute t-SNE
    tsne = TSNE(n_components=2, random_state=42)
    tsne_results = tsne.fit_transform(df.values)
    plot_projection(tsne_results,
                    markers=markers,
                    colors=colors,
                    save_path= save_path,
                    title=title,
                    name=name,
                    legend=legend)


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
            '2_way_norm': 'Two-way normalized',
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
name_map ={
    'robust': 'Robust',
    'standardized': 'Standardized',
    'robust+': 'Robust Norm\ntwo ways',
    'standardized+': 'Standardized\ntwo ways',
    # '2_way_norm_og': 'Study corrected\nthen genes Robust Normalized OG',
    '2_way_norm': 'Two-way\nNormalized',
    'study_corrected': 'Study\nCorrected',
    'imputed': 'No correction'
}

def plot_summary_scores_modified(scores_dict: dict, title: str, file_name: str, output_dir: str):
    """
    Generates an improved bar plot with grouped section labels for clustering scores,
    with bars re-ordered to group TREATMENT, TISSUE, and study scores consecutively
    for each matrix type. Includes a legend for TISSUE, TREATMENT, and study.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'bar_rand_ind_grouped.svg').
        output_dir (str): The directory to save the plot in.
    """
    # Note: This function assumes 're', 'os', and 'matplotlib.pyplot as plt' 
    # are imported in the global scope where this function is defined.
    
    # Import required for custom legend patches
    import matplotlib.patches as mpatches
    
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
        'study': '#d62728',      # Red (as per hex, though comment said Gray)
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

    # This assumes 'name_map' exists in the scope where the function is called
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
        ax.text(mid_point, -0.08, label_text, ha='center', va='top', fontsize=16, weight='bold',
                transform=ax.get_xaxis_transform(), rotation=0)

    # --- 5. Final Plot Adjustments ---
    
    # **MODIFICATION**: Remove the x-axis tick marks and labels
    ax.set_xticks([])
    
    # **MODIFICATION**: Add a custom legend
    treatment_patch = mpatches.Patch(color=color_map['treatment'], label='TREATMENT')
    tissue_patch = mpatches.Patch(color=color_map['tissue'], label='TISSUE')
    study_patch = mpatches.Patch(color=color_map['study'], label='study')
    ax.legend(handles=[treatment_patch, tissue_patch, study_patch], fontsize=12)

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
    Includes a legend for TISSUE and TREATMENT colors.

    Args:
        scores_dict (dict): Dictionary with labels as keys and scores as values.
        title (str): The title for the plot.
        file_name (str): The name of the output file (e.g., 'bar_rand_ind_relative.svg').
        output_dir (str): The directory to save the plot in.
    """
    # Note: This function assumes 're', 'os', and 'matplotlib.pyplot as plt' 
    # are imported in the global scope where this function is defined.
    
    # Import required for custom legend patches
    import matplotlib.patches as mpatches
    
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
    # This assumes 'name_map' exists in the scope where the function is called
    for start, end in zip(section_starts, section_ends):
        # Draw a vertical line at the start of each new section (except the first)
        if start > 0:
            ax.axvline(x=start - 0.5, color='gray', linestyle='--', linewidth=1.5)

        # Calculate the midpoint of the section to place the group label
        mid_point = start + (end - 1 - start) / 2
        
        raw_label_text = processing_types[start] 
        label_text = name_map.get(raw_label_text, raw_label_text)
        
        # Add the group label text below the x-axis
        ax.text(mid_point, -0.08, label_text, ha='center', va='top', fontsize=16, weight='bold',
                transform=ax.get_xaxis_transform(), rotation=0)

    # --- 5. Final Plot Adjustments ---
    
    # **MODIFICATION**: Remove the x-axis tick marks and labels
    ax.set_xticks([])
    
    # **MODIFICATION**: Add a custom legend for TISSUE and TREATMENT
    treatment_patch = mpatches.Patch(color=color_map['treatment'], label='TREATMENT')
    tissue_patch = mpatches.Patch(color=color_map['tissue'], label='TISSUE')
    ax.legend(handles=[treatment_patch, tissue_patch], fontsize=14)

    # **MODIFICATION**: Update Y-label (more robustly)
    ax.set_ylabel('Relative Score (X - Study Score)', fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Ensure the output directory exists and save the figure
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, file_name),bbox_inches='tight')
    plt.close()
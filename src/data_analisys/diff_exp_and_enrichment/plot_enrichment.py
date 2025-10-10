import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
from adjustText import adjust_text
from typing import Optional
import sys
module_dir = './'
sys.path.append(module_dir)
from src.constants import *

def plot_enrichment_results_old(results_df: pd.DataFrame,stress):
    """
    Generates a comparative dot plot for GSEA enrichment results.

    The plot visualizes the Normalized Enrichment Score (NES) and the
    statistical significance (FDR q-value) for multiple gene sets.

    Args:
        results_df (pd.DataFrame): DataFrame containing GSEA results with columns
                                   'Term', 'NES', and 'FDR q-val'.
    """
    # 1. Prepare the data for plotting
    df_plot = results_df[['Term', 'NES', 'FDR q-val']].copy()

    # Clean up the term names for better readability on the y-axis
    # Extracts "response to heat" from "GO:0009408 (response to heat)"
    df_plot['Gene Set'] = df_plot['Term'].str.extract(r'\((.*?)\)')[0].str.capitalize()

    # Calculate -log10(FDR) for dot size. This makes smaller FDR values (more significant)
    # result in larger dots. We add a small constant to handle FDR=0 cases.
    min_fdr = df_plot[df_plot['FDR q-val'] > 0]['FDR q-val'].min()
    df_plot['-log10(FDR)'] = -np.log10(np.array(df_plot['FDR q-val'].replace(0.0, min_fdr * 0.1),dtype=float))

    # Sort values by NES for a more organized plot
    df_plot = df_plot.sort_values(by='NES', ascending=True)

    # 2. Create the plot
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))

    # Use a diverging color palette to distinguish positive and negative enrichment
    palette = sns.color_palette("vlag_r", as_cmap=True)

    # Create the scatter plot, encoding NES as color and -log10(FDR) as size
    sns.scatterplot(
        data=df_plot,
        x='NES',
        y='Gene Set',
        size='-log10(FDR)',
        hue='NES',
        palette=palette,
        sizes=(50, 400),  # Min and max dot sizes
        legend='auto',
        ax=ax,
        edgecolor='black',
        linewidth=0.5
    )

    # 3. Customize the plot for clarity
    # Add a vertical line at x=0 to separate positive and negative enrichment
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=1)

    # Set titles and labels
    ax.set_title(f'GSEA Enrichment Results for {stress} - Stress Responses', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    ax.set_ylabel('')  # Y-axis label is not needed as terms are descriptive

    # Customize the legend
    handles, labels = ax.get_legend_handles_labels()
    # Separate size and color legends
    size_legend_index = df_plot.columns.get_loc('-log10(FDR)') + 1 # type: ignore
    color_legend_index = df_plot.columns.get_loc('NES') + 1 # type: ignore

    # Place legends neatly
    size_legend = ax.legend(handles[1:size_legend_index], labels[1:size_legend_index],
                            title='-log10(FDR)', bbox_to_anchor=(1.02, 0.7),
                            loc='upper left', labelspacing=1.2)
    ax.add_artist(size_legend)
    color_legend = ax.legend(handles[size_legend_index:], labels[size_legend_index:],
                             title='NES', bbox_to_anchor=(1.02, 0.4),
                             loc='upper left')
    ax.add_artist(color_legend)
    # Improve layout
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # type: ignore # Adjust rect to make space for legends
    plt.savefig(f'test/{stress}.svg')
    plt.close()

def plot_enrichment_scatter_all_labels(enrichment_df: pd.DataFrame, title: str = "Gene Set Enrichment Analysis", save_path: str = None):
    """
    Generates a scatter plot from gene set enrichment analysis results with non-overlapping labels.

    The plot visualizes:
    - X-axis: Normalized Enrichment Score (NES)
    - Y-axis: -log10 of the nominal p-value
    - Color: -log10 of the FDR q-value

    Args:
        enrichment_df (pd.DataFrame): A DataFrame containing the enrichment results.
                                      Must include the columns 'NES', 'NOM p-val', 'FDR q-val', and 'Term'.
        title (str, optional): The title for the plot.
        save_path (str, optional): File path to save the plot image.
    """
    # --- 1. Data Preparation ---
    df = enrichment_df.copy()
    df['-log10_pval'] = -np.log10(df['NOM p-val'].astype(float).replace(0, 1e-300))
    df['-log10_qval'] = -np.log10(df['FDR q-val'].astype(float).replace(0, 1e-300))

    # --- 2. Plotting ---
    fig, ax = plt.subplots(figsize=(12, 8))
    scatter = ax.scatter(
        x=df['NES'],
        y=df['-log10_pval'],
        c=df['-log10_qval'],
        cmap='viridis',
        s=50,
        alpha=0.8
    )

    # --- 3. Labels and Titles ---
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=12)
    ax.set_ylabel("-log10(Nominal P-value)", fontsize=12)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=1)
    
    significance_threshold = -np.log10(0.05)
    ax.axhline(y=significance_threshold, color='red', linestyle=':', linewidth=1.5, label='p = 0.05')
    ax.legend()

    # --- 4. Add labels with overlap prevention ---
    texts = []
    # To avoid clutter, we'll label the top 20 most significant terms (by FDR)
    df_to_label = df#.nsmallest(20, 'FDR q-val')

    for index, row in df_to_label.iterrows():
        term_label = row['Term'].split('(')[-1].replace(')', '').strip()
        texts.append(ax.text(
            x=row['NES'],
            y=row['-log10_pval'],
            s=term_label,
            fontdict={'size': 9}
        ))

    # Use adjust_text to automatically position labels to avoid overlap
    # TODO: only adjust the text that is at the pereto fronteer of Normalized enrichment score, -log10_pval p-value and -log10_qval
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    # --- 5. Color Bar ---
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("-log10(FDR q-value)", fontsize=12)

    # --- 6. Final Touches ---
    plt.tight_layout()

    if save_path:
        save_path = f'plots_enrichment/{save_path}'
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")
    else:
        plt.show()
    plt.close()


def plot_enrichment_results(results_df: pd.DataFrame,stress):
    """
    Generates a comparative dot plot for GSEA enrichment results.

    The plot visualizes the Normalized Enrichment Score (NES) and the
    statistical significance (FDR q-value) for multiple gene sets.

    Args:
        results_df (pd.DataFrame): DataFrame containing GSEA results with columns
                                   'Term', 'NES', and 'FDR q-val'.
    """
    # 1. Prepare the data for plotting
    df_plot = results_df[['Term', 'NES', 'FDR q-val']].copy()

    # Clean up the term names for better readability on the y-axis
    # Extracts "response to heat" from "GO:0009408 (response to heat)"
    df_plot['Gene Set'] = df_plot['Term'].str.extract(r'\((.*?)\)')[0].str.capitalize()

    # Calculate -log10(FDR) for dot size. This makes smaller FDR values (more significant)
    # result in larger dots. We add a small constant to handle FDR=0 cases.
    min_fdr = df_plot[df_plot['FDR q-val'] > 0]['FDR q-val'].min()
    df_plot['-log10(FDR)'] = -np.log10(np.array(df_plot['FDR q-val'].replace(0.0, min_fdr * 0.1),dtype=float))

    # Sort values by NES for a more organized plot
    df_plot = df_plot.sort_values(by='NES', ascending=True)

    # 2. Create the plot
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax = plt.subplots(figsize=(10, 7))

    # Use a diverging color palette to distinguish positive and negative enrichment
    palette = sns.color_palette("vlag_r", as_cmap=True)

    # Create the scatter plot, encoding NES as color and -log10(FDR) as size
    sns.scatterplot(
        data=df_plot,
        x='NES',
        y='Gene Set',
        size='-log10(FDR)',
        hue='NES',
        palette=palette,
        sizes=(50, 400),  # Min and max dot sizes
        legend='auto',
        ax=ax,
        edgecolor='black',
        linewidth=0.5
    )

    # 3. Customize the plot for clarity
    # Add a vertical line at x=0 to separate positive and negative enrichment
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=1)

    # Set titles and labels
    ax.set_title(f'GSEA Enrichment Results for {stress} - Stress Responses', fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=12)
    ax.set_ylabel('')  # Y-axis label is not needed as terms are descriptive

    # Customize the legend
    handles, labels = ax.get_legend_handles_labels()
    # Separate size and color legends
    size_legend_index = df_plot.columns.get_loc('-log10(FDR)') + 1
    color_legend_index = df_plot.columns.get_loc('NES') + 1

    # Place legends neatly
    size_legend = ax.legend(handles[1:size_legend_index], labels[1:size_legend_index],
                            title='-log10(FDR)', bbox_to_anchor=(1.02, 0.7),
                            loc='upper left', labelspacing=1.2)
    ax.add_artist(size_legend)
    color_legend = ax.legend(handles[size_legend_index:], labels[size_legend_index:],
                             title='NES', bbox_to_anchor=(1.02, 0.4),
                             loc='upper left')
    ax.add_artist(color_legend)
    # Improve layout
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust rect to make space for legends
    plt.savefig(f'test/{stress}.svg')
    plt.close()

def find_pareto_frontier_indices(df: pd.DataFrame, margin: float = 0.05) -> pd.Index:
    """
    Identifies indices of points on or near the Pareto frontier.

    A point is on the frontier if it's not dominated by any other point.
    This version includes a margin, so a point is only considered dominated if
    another point is significantly better (i.e., outside the margin).

    We aim to maximize: abs(NES), -log10_pval, -log10_qval.

    Args:
        df (pd.DataFrame): DataFrame with columns 'NES', '-log10_pval', '-log10_qval'.
        margin (float): The relative margin (e.g., 0.05 for 5%) to consider points
                        as being "close" to the frontier.

    Returns:
        pd.Index: Indices of rows on or near the Pareto frontier.
    """
    values = df[['NES', '-log10_pval', '-log10_qval']].copy()
    values['NES'] = values['NES'].abs()
    points = values.to_numpy()
    
    num_points = points.shape[0]
    is_on_frontier = np.ones(num_points, dtype=bool)

    for i in range(num_points):
        for j in range(num_points):
            if i == j:
                continue
            
            # A point j dominates point i if it is better than i * (1 + margin) in all criteria
            # and strictly better in at least one. By adding the margin to point i, we make it
            # "harder" for j to dominate it, thus keeping points close to the frontier.
            if np.all(points[j] >= points[i] * (1 + margin)) and np.any(points[j] > points[i]):
                is_on_frontier[i] = False
                break
    
    return df.index[is_on_frontier]

def plot_enrichment_scatter(enrichment_df: pd.DataFrame, title: str = "Gene Set Enrichment Analysis", save_path: Optional[str] =None):
    """
    Generates a scatter plot from gene set enrichment analysis results with non-overlapping labels.

    The plot visualizes:
    - X-axis: Normalized Enrichment Score (NES)
    - Y-axis: -log10 of the nominal p-value
    - Color: -log10 of the FDR q-value

    Args:
        enrichment_df (pd.DataFrame): A DataFrame containing the enrichment results.
                                      Must include the columns 'NES', 'NOM p-val', 'FDR q-val', and 'Term'.
        title (str, optional): The title for the plot.
        save_path (str, optional): File path to save the plot image.
    """
    # --- 1. Data Preparation ---
    df = enrichment_df.copy()
    df['-log10_pval'] = -np.log10(df['NOM p-val'].astype(float).replace(0, 1e-3))
    df['-log10_qval'] = -np.log10(df['FDR q-val'].astype(float).replace(0, 1e-3))

    # --- 2. Plotting ---
    fig, ax = plt.subplots(figsize=(12, 8))
    scatter = ax.scatter(
        x=df['NES'],
        y=df['-log10_pval'],
        c=df['-log10_qval'],
        cmap='viridis',
        s=50,
        alpha=0.8
    )

    # --- 3. Labels and Titles ---
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=12)
    ax.set_ylabel("-log10(Nominal P-value)", fontsize=12)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.axvline(x=0, color='grey', linestyle='--', linewidth=1)
    
    significance_threshold = -np.log10(0.05)
    ax.axhline(y=significance_threshold, color='red', linestyle=':', linewidth=1.5, label='p = 0.05')
    ax.legend()

    # --- 4. Add labels with overlap prevention ---
    texts = []
    
    # Find points on or near the Pareto frontier to label.
    # The margin parameter (e.g., 0.05) includes points within 5% of the frontier.
    pareto_indices = find_pareto_frontier_indices(df, margin=0.01)
    df_to_label = df.loc[pareto_indices]

    for index, row in df_to_label.iterrows():
        term_label = row['Term'].split('(')[-1].replace(')', '').strip()
        texts.append(ax.text(
            x=row['NES'],
            y=row['-log10_pval'],
            s=term_label,
            fontdict={'size': 9}
        ))

    # Use adjust_text to automatically position labels to avoid overlap
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    # --- 5. Color Bar ---
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("-log10(FDR q-value)", fontsize=12)

    # --- 6. Final Touches ---
    plt.tight_layout()

    if save_path:
        save_path = f'plots_enrichment/{save_path}'
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to {save_path}")
    else:
        plt.show()
    plt.close()

###################
import plotly.graph_objects as go
import uuid
import json
def plot_enrichment_scatter_interactive(enrichment_df: pd.DataFrame, title: str = "Gene Set Enrichment Analysis", save_path: str = "interactive_plot.html"):
    """
    Generates a self-contained, interactive HTML scatter plot from GSEA results.

    This version embeds custom JavaScript to create a searchable, multi-select dropdown
    menu within a single, portable HTML file.

    Features:
    - A searchable dropdown to find and highlight multiple gene sets.
    - A button to toggle labels for the most significant terms.
    - Zoom, pan, and hover over points for detailed information.

    Args:
        enrichment_df (pd.DataFrame): DataFrame with GSEA results.
        title (str, optional): The title for the plot.
        save_path (str, optional): The path to save the output HTML file.
    """
    # --- 1. Data Preparation ---
    df = enrichment_df.copy()
    df['-log10_pval'] = -np.log10(df['NOM p-val'].astype(float).replace(0, 1e-10))
    df['-log10_qval'] = -np.log10(df['FDR q-val'].astype(float).replace(0, 1e-10))
    
    df['hover_text'] = df.apply(
        lambda row: f"""<b>{row['Term']}</b><br><br>
NES: {row['NES']:.3f}<br>
Nominal p-val: {row['NOM p-val']:.3g}<br>
FDR q-val: {row['FDR q-val']:.3g}
""",
        axis=1
    )
    
    # Identify significant terms to label
    pareto_indices = find_pareto_frontier_indices(df)
    df_to_label = df.loc[pareto_indices]

    # --- 2. Create the Plotly Figure ---
    fig = go.Figure()
    plot_div_id = f'plotly-graph-{uuid.uuid4()}' # Unique ID for the graph div

    # Trace 0: Main scatter plot
    fig.add_trace(go.Scatter(
        x=df['NES'], y=df['-log10_pval'], mode='markers', hoverinfo='text',
        hovertext=df['hover_text'], name='Gene Sets',
        marker=dict(
            color=df['-log10_qval'], colorscale='Viridis', showscale=True,
            colorbar=dict(title="-log10(FDR)"), size=8, symbol='circle'
        )
    ))

    # Trace 1: Labels for significant terms
    fig.add_trace(go.Scatter(
        x=df_to_label['NES'], y=df_to_label['-log10_pval'], mode='text',
        text=df_to_label['Term'], textposition="top right",
        textfont=dict(size=10, color='#444'), hoverinfo='none', name='Labels',
        visible=True # Initially visible
    ))

    # Trace 2: Placeholder for highlighted points (will be controlled by JavaScript)
    fig.add_trace(go.Scatter(
        x=[], y=[], mode='markers', hoverinfo='none', name='Selected', showlegend=False,
        marker=dict(color='red', size=16, symbol='star', line=dict(width=1, color='black'))
    ))
    
    # --- 3. Configure Layout and Controls ---
    fig.update_layout(
        title=dict(text=f'<b>{title}</b>', x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title="-log10(Nominal P-value)",
        template='plotly_white', height=700, hovermode='closest',
        updatemenus=[
            dict(
                type="buttons", direction="right", active=0, x=0.9, y=1.1,
                buttons=[
                    dict(label="Show Labels", method="restyle", args=[{"visible": [True, True, True]}, [0, 1, 2]]),
                    dict(label="Hide Labels", method="restyle", args=[{"visible": [True, False, True]}, [0, 1, 2]]),
                ]
            )
        ]
    )
    fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="grey")
    fig.add_hline(y=-np.log10(0.05), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text="p = 0.05", annotation_position="bottom right")
    fig.add_hline(y=-np.log10(1/2000), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text=f"p = {1/2000}", annotation_position="bottom right")

    # --- 4. Prepare Data and HTML for Injection ---
    
    # Generate the core Plotly graph HTML
    plot_div = fig.to_html(
        full_html=False, 
        include_plotlyjs='cdn', 
        div_id=plot_div_id
    )

    # Prepare data for JavaScript: a map of Term -> Coordinates
    coords_data = df[['Term', 'NES', '-log10_pval']].to_dict(orient='records')
    coords_json = json.dumps(coords_data)
    
    # Create HTML <option> elements for the dropdown
    options_html = "".join([f'<option value="{term}">{term}</option>' for term in sorted(df['Term'])])

    # --- 5. Construct the Final HTML Document ---
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8" />
        <title>{title}</title>
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/choices.js/public/assets/styles/choices.min.css"/>
        <script src="https://cdn.jsdelivr.net/npm/choices.js/public/assets/scripts/choices.min.js"></script>
        <style>
            body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }}
            .container {{ max-width: 1200px; margin: 20px auto; text-align: center; }}
            .choices {{ margin: 0 auto 20px auto; max-width: 80%; text-align: left; }}
            .choices__inner {{ background-color: #f9f9f9; }}
        </style>
    </head>
    <body>
        <div class="container">
            {plot_div}
            
            <label for="term-select">Search and select gene sets to highlight:</label>
            <select id="term-select" multiple>
                {options_html}
            </select>
        </div>

        <script>
            // --- JavaScript to link the dropdown with the Plotly graph ---
            
            // 1. Data passed from Python
            const termData = {coords_json};
            const plotDivId = '{plot_div_id}';

            // 2. Create a Map for fast coordinate lookups (Term -> {{NES, -log10_pval}})
            const coordMap = new Map(termData.map(item => [item.Term, item]));

            // 3. Initialize the searchable dropdown using Choices.js
            const choices = new Choices('#term-select', {{
                removeItemButton: true,
                searchResultLimit: 150,
                shouldSort: false, // Options are pre-sorted in Python
                placeholder: true,
                placeholderValue: 'Type to search...',
            }});

            // 4. Listen for changes in the dropdown selection
            document.getElementById('term-select').addEventListener('change', function(event) {{
                const selectedTerms = choices.getValue(true); // Get selected values as an array
                
                const highlightX = [];
                const highlightY = [];

                // Find coordinates for each selected term
                selectedTerms.forEach(term => {{
                    const data = coordMap.get(term);
                    if (data) {{
                        highlightX.push(data.NES);
                        highlightY.push(data['-log10_pval']);
                    }}
                }});

                // 5. Update the 'Selected' trace (index 2) on the Plotly graph
                Plotly.restyle(plotDivId, {{
                    x: [highlightX],
                    y: [highlightY]
                }}, [2]);
            }});

            // Move the dropdown above the plot for better layout
            const plotContainer = document.getElementById(plotDivId).parentElement;
            const dropdown = document.querySelector('.choices').parentElement;
            plotContainer.parentNode.insertBefore(dropdown, plotContainer);
            
        </script>
    </body>
    </html>
    """

    # --- 6. Save to HTML file ---
    if save_path:
        output_dir = os.path.dirname(save_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write(html_template)
        print(f"Interactive plot saved to: {os.path.abspath(save_path)}")

############
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
import json
from typing import Optional, Dict, List, Any

def find_pareto_frontier_indices(df: pd.DataFrame, margin: float = 0.0) -> pd.Index:
    """
    Identifies indices of points on or near the Pareto frontier.
    We aim to maximize: abs(NES), -log10_pval, -log10_qval.
    """
    values = df[['NES', '-log10_pval', '-log10_qval']].copy()
    values['NES'] = values['NES'].abs()
    points = values.to_numpy()
    
    num_points = points.shape[0]
    is_on_frontier = np.ones(num_points, dtype=bool)

    for i in range(num_points):
        for j in range(num_points):
            if i == j:
                continue
            
            if np.all(points[j] >= points[i] * (1 + margin)) and np.any(points[j] > points[i]):
                is_on_frontier[i] = False
                break
    
    return df.index[is_on_frontier]

def aggregate_term_trajectories(base_results_dir: str) -> Dict[str, List[Dict[str, Any]]]:
    """
    Scans the directory structure for all GSEA result files and aggregates
    the data to track each GO term's performance across all conditions.
    
    Args:
        base_results_dir: The root directory where all GSEA results are stored.
                          e.g., 'GSEA_results/'
    
    Returns:
        A dictionary where keys are GO Terms and values are a list of their
        results across all found experimental conditions.
    """
    trajectory_data = {}
    print(f"Scanning for GSEA results in '{base_results_dir}'...")

    for root, _, files in os.walk(base_results_dir):
        for file in files:
            # Assuming the raw results are stored as CSV files.
            if file.endswith('.csv'):
                file_path = os.path.join(root, file)
                try:
                    # Extract parameters from the path relative to the base directory
                    relative_path = os.path.relpath(file_path, base_results_dir)
                    path_parts = relative_path.split(os.sep)
                    
                    if len(path_parts) < 7: continue # Skip files not in the expected structure

                    # Example structure: 1.0/full/All-Tissues/corrected/0/pure/Cold Stress.csv
                    params = {
                        'version': path_parts[0],
                        'dataset': path_parts[1],
                        'tissue': path_parts[2],
                        'normalization': path_parts[3],
                        'param_value': path_parts[4],
                        'purity': path_parts[5],
                        'stress': path_parts[6].replace('.csv', '')
                    }

                    df = pd.read_csv(file_path)
                    for _, row in df.iterrows():
                        term = row['Term']
                        if term not in trajectory_data:
                            trajectory_data[term] = []
                        
                        result_entry = params.copy()
                        result_entry['NES'] = row['NES']
                        result_entry['NOM p-val'] = row['NOM p-val']
                        result_entry['FDR q-val'] = row['FDR q-val']
                        trajectory_data[term].append(result_entry)

                except Exception as e:
                    print(f"Warning: Could not process file '{file_path}'. Error: {e}")

    print(f"Finished scanning. Found trajectory data for {len(trajectory_data)} unique terms.")
    return trajectory_data


def build_interactive_explorer(
    current_enrichment_df: pd.DataFrame, 
    all_trajectories_data: dict,
    title: str = "Gene Set Enrichment Analysis", 
    save_path: str = "interactive_plot.html"
):
    """
    Generates a self-contained, interactive HTML file for exploring GSEA results.
    """
    # --- 1. Prepare data for the main plot ---
    df = current_enrichment_df.copy()
    df['-log10_pval'] = -np.log10(df['NOM p-val'].astype(float).replace(0, 1e-300))
    df['-log10_qval'] = -np.log10(df['FDR q-val'].astype(float).replace(0, 1e-300))
    df['hover_text'] = df.apply(
        lambda row: f"<b>{row['Term']}</b><br>Click to select and track this term", axis=1
    )

    # --- 2. Serialize all trajectory data for embedding ---
    # This will be a large JSON string embedded in the HTML
    trajectory_json = json.dumps(all_trajectories_data)

    # --- 3. Generate the main Plotly figure ---
    fig = go.Figure()

    # Main scatter plot trace
    fig.add_trace(go.Scatter(
        x=df['NES'], y=df['-log10_pval'], mode='markers',
        hoverinfo='text', hovertext=df['hover_text'],
        customdata=df['Term'], # Store Term name in customdata to access on click
        marker=dict(
            color=df['-log10_qval'], colorscale='Viridis', showscale=True,
            colorbar=dict(title="-log10(FDR q-val)"), size=8
        )
    ))

    # Trace for the highlighted point
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode='markers',
        marker=dict(color='red', size=16, symbol='star', line=dict(width=1, color='black')),
        hoverinfo='none', name='Selected'
    ))

    # --- 4. HTML and JavaScript Generation ---
    # This is the core of the new functionality
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{title}</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ font-family: Arial, sans-serif; display: flex; }}
            #main-plot {{ width: 60%; }}
            #trajectory-dashboard {{ width: 40%; padding-left: 20px; }}
            .trajectory-plot {{ height: 250px; width: 100%; }}
            h2, h3 {{ color: #333; }}
        </style>
    </head>
    <body>
        <div id="main-plot"></div>
        <div id="trajectory-dashboard">
            <h2>Term Trajectory Explorer</h2>
            <h3 id="selected-term-title">Select a point on the plot to begin.</h3>
            <div id="trajectory-plot-nes" class="trajectory-plot"></div>
            <div id="trajectory-plot-fdr" class="trajectory-plot"></div>
        </div>

        <script>
            // Embed the main plot and all trajectory data
            const mainPlotData = {fig.to_json()};
            const trajectoryData = JSON.parse('{trajectory_json}');
            let selectedTerm = null;

            // Render the main plot
            Plotly.newPlot('main-plot', mainPlotData);
            const mainPlotDiv = document.getElementById('main-plot');

            // --- Event Listener for Clicks ---
            mainPlotDiv.on('plotly_click', function(data){{
                if(data.points.length > 0) {{
                    selectedTerm = data.points[0].customdata;
                    document.getElementById('selected-term-title').innerText = "Trajectory for: " + selectedTerm;
                    
                    // Highlight the selected point
                    const update = {{
                        x: [[data.points[0].x]],
                        y: [[data.points[0].y]]
                    }};
                    Plotly.restyle(mainPlotDiv, update, [1]); // Update the highlight trace

                    // Draw the trajectory plots
                    drawTrajectoryPlots(selectedTerm);
                }}
            }});

            // --- Function to Draw Trajectory Plots ---
            function drawTrajectoryPlots(term) {{
                const termData = trajectoryData[term];
                if (!termData) return;

                // --- Plot 1: NES vs. Parameter Value (0, 10, 15) ---
                const nesTrace = {{
                    x: termData.map(d => d.param_value),
                    y: termData.map(d => d.NES),
                    mode: 'markers',
                    type: 'scatter',
                    text: termData.map(d => `Stress: ${{d.stress}}<br>Tissue: ${{d.tissue}}`),
                    marker: {{ size: 10 }}
                }};
                const nesLayout = {{
                    title: 'NES vs. Parameter Value',
                    xaxis: {{ title: 'Parameter (0, 10, 15)' }},
                    yaxis: {{ title: 'Normalized Enrichment Score (NES)' }}
                }};
                Plotly.newPlot('trajectory-plot-nes', [nesTrace], nesLayout);

                // --- Plot 2: -log10(FDR) vs. Stress Type ---
                const fdrTrace = {{
                    x: termData.map(d => d.stress),
                    y: termData.map(d => -Math.log10(d['FDR q-val'] || 1e-10)),
                    mode: 'markers',
                    type: 'scatter',
                    text: termData.map(d => `Tissue: ${{d.tissue}}<br>Purity: ${{d.purity}}`),
                    marker: {{ size: 10 }}
                }};
                const fdrLayout = {{
                    title: '-log10(FDR) vs. Stress Type',
                    xaxis: {{ title: 'Stress Type', tickangle: -45 }},
                    yaxis: {{ title: '-log10(FDR q-val)' }}
                }};
                Plotly.newPlot('trajectory-plot-fdr', [fdrTrace], fdrLayout);
            }}
        </script>
    </body>
    </html>
    """

    # --- 5. Save the HTML file ---
    if save_path:
        output_dir = os.path.dirname(save_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write(html_template)
        print(f"Interactive explorer saved to: {save_path}")


##### NEW plot method:
def plot_enrichment_scatter_interactive_2(enrichment_df: pd.DataFrame, title: str = "Gene Set Enrichment Analysis", save_path: str = "interactive_plot.html",treatments=None):
    """
    Generates a self-contained, interactive HTML scatter plot from GSEA results.

    This version embeds custom JavaScript to create a searchable, multi-select dropdown
    menu within a single, portable HTML file, plus navigation buttons that change URLs
    while preserving highlighted terms.

    Features:
    - A searchable dropdown to find and highlight multiple gene sets.
    - A button to toggle labels for the most significant terms.
    - Zoom, pan, and hover over points for detailed information.
    - Navigation buttons that change URLs while preserving highlighted terms.

    Args:
        enrichment_df (pd.DataFrame): DataFrame with GSEA results.
        title (str, optional): The title for the plot.
        save_path (str, optional): The path to save the output HTML file.
    """
    # --- 1. Data Preparation ---
    df = enrichment_df.copy()
    df['-log10_pval'] = -np.log10(df['NOM p-val'].astype(float).replace(0, 1e-10))
    df['-log10_qval'] = -np.log10(df['FDR q-val'].astype(float).replace(0, 1e-10))
    
    df['hover_text'] = df.apply(
        lambda row: f"""<b>{row['Term']}</b><br><br>
NES: {row['NES']:.3f}<br>
Nominal p-val: {row['NOM p-val']:.3g}<br>
FDR q-val: {row['FDR q-val']:.3g}
""",
        axis=1
    )
    
    # Identify significant terms to label
    pareto_indices = find_pareto_frontier_indices(df)
    df_to_label = df.loc[pareto_indices]

    # --- 2. Create the Plotly Figure ---
    fig = go.Figure()
    plot_div_id = f'plotly-graph-{uuid.uuid4()}' # Unique ID for the graph div

    # Trace 0: Main scatter plot
    fig.add_trace(go.Scatter(
        x=df['NES'], y=df['-log10_pval'], mode='markers', hoverinfo='text',
        hovertext=df['hover_text'], name='Gene Sets',
        marker=dict(
            color=df['-log10_qval'], colorscale='Viridis', showscale=True,
            colorbar=dict(title="-log10(FDR)"), size=8, symbol='circle'
        )
    ))

    # Trace 1: Labels for significant terms
    fig.add_trace(go.Scatter(
        x=df_to_label['NES'], y=df_to_label['-log10_pval'], mode='text',
        text=df_to_label['Term'], textposition="top right",
        textfont=dict(size=10, color='#444'), hoverinfo='none', name='Labels',
        visible=True # Initially visible
    ))

    # Trace 2: Placeholder for highlighted points (will be controlled by JavaScript)
    fig.add_trace(go.Scatter(
        x=[], y=[], mode='markers', hoverinfo='none', name='Selected', showlegend=False,
        marker=dict(color='red', size=16, symbol='star', line=dict(width=1, color='black'))
    ))
    
    # --- 3. Configure Layout and Controls ---
    fig.update_layout(
        title=dict(text=f'<b>{title}</b>', x=0.5),
        xaxis_title="Normalized Enrichment Score (NES)",
        yaxis_title="-log10(Nominal P-value)",
        template='plotly_white', height=700, hovermode='closest',
        updatemenus=[
            dict(
                type="buttons", direction="right", active=0, x=0.9, y=1.1,
                buttons=[
                    dict(label="Show Labels", method="restyle", args=[{"visible": [True, True, True]}, [0, 1, 2]]),
                    dict(label="Hide Labels", method="restyle", args=[{"visible": [True, False, True]}, [0, 1, 2]]),
                ]
            )
        ]
    )
    fig.add_vline(x=0, line_width=1, line_dash="dash", line_color="grey")
    fig.add_hline(y=-np.log10(0.05), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text="p = 0.05", annotation_position="bottom right")
    fig.add_hline(y=-np.log10(1/2000), line_width=1.5, line_dash="dot", line_color="red",
                  annotation_text=f"p = {1/2000}", annotation_position="bottom right")

    # --- 4. Prepare Data and HTML for Injection ---
    
    # Generate the core Plotly graph HTML
    plot_div = fig.to_html(
        full_html=False, 
        include_plotlyjs='cdn', 
        div_id=plot_div_id
    )

    # Prepare data for JavaScript: a map of Term -> Coordinates
    coords_data = df[['Term', 'NES', '-log10_pval']].to_dict(orient='records')
    coords_json = json.dumps(coords_data)
    
    # Create HTML <option> elements for the dropdown
    options_html = "".join([f'<option value="{term}">{term}</option>' for term in sorted(df['Term'])])

    # --- 5. Extract path components for navigation ---
    path_parts = save_path.split('/')
    
    # Define possible values for each parameter (you can customize these)
    dataset_types = ['full', 'sanity']
    tissues = ['All-Tissues', 'leaf']  # Add your actual tissues
    normalizations = ['2_way_norm', 'imputed', 'study_corrected']
    thresholds = ['0', '10', '15']
    purities = ['pure', 'mixed']
    stress_names = treatments  # Add your actual stress names
    
    # Initialize with default values
    version = '0.0'
    dataset_type = 'full'
    tissue = 'All-Tissues'
    normalization = '2_way_norm'
    threshold = '14'
    purity = 'pure'
    stress_name = 'Heat'
    
    # Try to extract values from path
    try:
        plots_idx = 2
        if len(path_parts) > plots_idx + 1:
            version = path_parts[plots_idx + 1]
        if len(path_parts) > plots_idx + 2:
            dataset_type = path_parts[plots_idx + 2]
        if len(path_parts) > plots_idx + 3:
            tissue = path_parts[plots_idx + 3]
        if len(path_parts) > plots_idx + 4:
            normalization = path_parts[plots_idx + 4]
        if len(path_parts) > plots_idx + 5:
            threshold = path_parts[plots_idx + 5]
        if len(path_parts) > plots_idx + 6:
            purity = path_parts[plots_idx + 6]
            
        # Extract stress name from filename
        if path_parts[-1].endswith(' Stress.html'):
            stress_name = path_parts[-1].replace('.html', '')
    except (ValueError, IndexError):
        pass  # Keep default values if extraction fails

    # Generate navigation URLs for each category
    nav_urls = {}
    
    # Helper function to generate URL
    def generate_url(v, dt, t, norm, thresh, p, stress):
        enrich_out = save_path.split('/')[2]
        return f"{GLOBAL_DIR_PATH}{FIGURES_DIR.split('.')[1][1:]}{enrich_out}/{v}/{dt}/{t}/{norm}/{thresh}/{p}/{stress}.html"
    
    # Generate URLs for each parameter type
    nav_urls['dataset_type'] = [
        (dt, generate_url(version, dt, tissue, normalization, threshold, purity, stress_name))
        for dt in dataset_types if dt != dataset_type
    ]
    
    nav_urls['tissue'] = [
        (t, generate_url(version, dataset_type, t, normalization, threshold, purity, stress_name))
        for t in tissues if t != tissue
    ]
    
    nav_urls['normalization'] = [
        (norm, generate_url(version, dataset_type, tissue, norm, threshold, purity, stress_name))
        for norm in normalizations if norm != normalization
    ]
    
    nav_urls['threshold'] = [
        (thresh, generate_url(version, dataset_type, tissue, normalization, thresh, purity, stress_name))
        for thresh in thresholds if thresh != threshold
    ]
    
    nav_urls['purity'] = [
        (p, generate_url(version, dataset_type, tissue, normalization, threshold, p, stress_name))
        for p in purities if p != purity
    ]
    
    nav_urls['stress'] = [
        (stress, generate_url(version, dataset_type, tissue, normalization, threshold, purity, stress))
        for stress in stress_names if stress != stress_name
    ]

    # --- 6. Construct the Final HTML Document ---
    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8" />
        <title>{title}</title>
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/choices.js/public/assets/styles/choices.min.css"/>
        <script src="https://cdn.jsdelivr.net/npm/choices.js/public/assets/scripts/choices.min.js"></script>
        <style>
            body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; }}
            .container {{ max-width: 1200px; margin: 20px auto; text-align: center; }}
            .choices {{ margin: 0 auto 20px auto; max-width: 80%; text-align: left; }}
            .choices__inner {{ background-color: #f9f9f9; }}
            .nav-buttons {{ margin: 10px 0; display: flex; justify-content: center; gap: 8px; flex-wrap: wrap; }}
            .nav-button {{ padding: 6px 12px; background-color: #4CAF50; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 12px; }}
            .nav-button:hover {{ background-color: #45a049; }}
            .nav-section {{ margin: 15px 0; padding: 10px; border: 1px solid #ddd; border-radius: 5px; }}
            .nav-section h3 {{ margin-bottom: 10px; color: #333; font-size: 14px; }}
            .current-info {{ background-color: #f0f8ff; padding: 10px; border-radius: 5px; margin: 15px 0; font-size: 14px; }}
            .no-buttons {{ color: #666; font-style: italic; font-size: 12px; }}
        </style>
    </head>
    <body>
        <div class="container">
            <!-- Current Plot Info -->
            <div class="current-info">
                <strong>Current Plot:</strong><br>
                Version: {version} | Dataset: {dataset_type} | Tissue: {tissue}<br>
                Normalization: {normalization} | Threshold: {threshold} | Purity: {purity} | Stress: {stress_name}
            </div>

            <!-- Navigation Sections -->
            <div class="nav-section">
                <h3>Navigate by Dataset Type</h3>
                <div class="nav-buttons" id="dataset-buttons">
                    {f'<span class="no-buttons">No other dataset types available</span>' if not nav_urls['dataset_type'] else ''}
                </div>
            </div>

            <div class="nav-section">
                <h3>Navigate by Tissue</h3>
                <div class="nav-buttons" id="tissue-buttons">
                    {f'<span class="no-buttons">No other tissues available</span>' if not nav_urls['tissue'] else ''}
                </div>
            </div>

            <div class="nav-section">
                <h3>Navigate by Normalization</h3>
                <div class="nav-buttons" id="normalization-buttons">
                    {f'<span class="no-buttons">No other normalizations available</span>' if not nav_urls['normalization'] else ''}
                </div>
            </div>

            <div class="nav-section">
                <h3>Navigate by Threshold</h3>
                <div class="nav-buttons" id="threshold-buttons">
                    {f'<span class="no-buttons">No other thresholds available</span>' if not nav_urls['threshold'] else ''}
                </div>
            </div>

            <div class="nav-section">
                <h3>Navigate by Purity</h3>
                <div class="nav-buttons" id="purity-buttons">
                    {f'<span class="no-buttons">No other purity types available</span>' if not nav_urls['purity'] else ''}
                </div>
            </div>

            <div class="nav-section">
                <h3>Navigate by Stress Type</h3>
                <div class="nav-buttons" id="stress-buttons">
                    {f'<span class="no-buttons">No other stress types available</span>' if not nav_urls['stress'] else ''}
                </div>
            </div>

            <!-- Main Plot -->
            {plot_div}
            
            <label for="term-select">Search and select gene sets to highlight:</label>
            <select id="term-select" multiple>
                {options_html}
            </select>
        </div>

        <script>
            // --- JavaScript to link the dropdown with the Plotly graph ---
            
            // 1. Data passed from Python
            const termData = {coords_json};
            const plotDivId = '{plot_div_id}';
            const navUrls = {json.dumps(nav_urls)};

            // 2. Create a Map for fast coordinate lookups (Term -> {{NES, -log10_pval}})
            const coordMap = new Map(termData.map(item => [item.Term, item]));

            // 3. Initialize the searchable dropdown using Choices.js
            const choices = new Choices('#term-select', {{
                removeItemButton: true,
                searchResultLimit: 150,
                shouldSort: false, // Options are pre-sorted in Python
                placeholder: true,
                placeholderValue: 'Type to search...',
            }});

            // 4. Function to create navigation buttons
            function createNavButtons(buttonData, containerId) {{
                const container = document.getElementById(containerId);
                // Clear any "no buttons" message
                if (container.querySelector('.no-buttons')) {{
                    container.innerHTML = '';
                }}
                
                buttonData.forEach(([label, url]) => {{
                    const button = document.createElement('button');
                    button.className = 'nav-button';
                    button.textContent = label;
                    button.onclick = function() {{
                        // Get currently selected terms
                        const selectedTerms = choices.getValue(true);
                        
                        // Navigate to new URL with selected terms as URL parameters
                        const params = new URLSearchParams();
                        selectedTerms.forEach(term => {{
                            params.append('highlight', term);
                        }});
                        
                        const separator = url.includes('?') ? '&' : '?';
                        window.location.href = url + separator + params.toString();
                    }};
                    
                    container.appendChild(button);
                }});
            }}

            // 5. Create all navigation buttons
            if (navUrls.dataset_type && navUrls.dataset_type.length > 0) {{
                createNavButtons(navUrls.dataset_type, 'dataset-buttons');
            }}

            if (navUrls.tissue && navUrls.tissue.length > 0) {{
                createNavButtons(navUrls.tissue, 'tissue-buttons');
            }}

            if (navUrls.normalization && navUrls.normalization.length > 0) {{
                createNavButtons(navUrls.normalization, 'normalization-buttons');
            }}

            if (navUrls.threshold && navUrls.threshold.length > 0) {{
                createNavButtons(navUrls.threshold, 'threshold-buttons');
            }}

            if (navUrls.purity && navUrls.purity.length > 0) {{
                createNavButtons(navUrls.purity, 'purity-buttons');
            }}

            if (navUrls.stress && navUrls.stress.length > 0) {{
                createNavButtons(navUrls.stress, 'stress-buttons');
            }}

            // 6. Listen for changes in the dropdown selection
            document.getElementById('term-select').addEventListener('change', function(event) {{
                const selectedTerms = choices.getValue(true); // Get selected values as an array
                
                const highlightX = [];
                const highlightY = [];

                // Find coordinates for each selected term
                selectedTerms.forEach(term => {{
                    const data = coordMap.get(term);
                    if (data) {{
                        highlightX.push(data.NES);
                        highlightY.push(data['-log10_pval']);
                    }}
                }});

                // Update the 'Selected' trace (index 2) on the Plotly graph
                Plotly.restyle(plotDivId, {{
                    x: [highlightX],
                    y: [highlightY]
                }}, [2]);

                // Store selected terms in sessionStorage to preserve across navigation
                sessionStorage.setItem('selectedTerms', JSON.stringify(selectedTerms));
            }});

            // 7. Restore selected terms from URL parameters or sessionStorage on page load
            function restoreSelectedTerms() {{
                const urlParams = new URLSearchParams(window.location.search);
                const highlightTerms = urlParams.getAll('highlight');
                
                if (highlightTerms.length > 0) {{
                    // Set selections from URL parameters
                    highlightTerms.forEach(term => {{
                        choices.setChoiceByValue(term);
                    }});
                    // Trigger the change event to update the plot
                    document.getElementById('term-select').dispatchEvent(new Event('change'));
                }} else {{
                    // Try to restore from sessionStorage
                    const storedTerms = sessionStorage.getItem('selectedTerms');
                    if (storedTerms) {{
                        const terms = JSON.parse(storedTerms);
                        terms.forEach(term => {{
                            choices.setChoiceByValue(term);
                        }});
                        document.getElementById('term-select').dispatchEvent(new Event('change'));
                    }}
                }}
            }}

            // Call restore function when page loads
            document.addEventListener('DOMContentLoaded', restoreSelectedTerms);

            // Move the dropdown above the plot for better layout
            const plotContainer = document.getElementById(plotDivId).parentElement;
            const dropdown = document.querySelector('.choices').parentElement;
            plotContainer.parentNode.insertBefore(dropdown, plotContainer);
            
        </script>
    </body>
    </html>
    """

    # --- 7. Save to HTML file ---
    if save_path:
        output_dir = os.path.dirname(save_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write(html_template)
        print(f"Interactive plot saved to: {os.path.abspath(save_path)}")
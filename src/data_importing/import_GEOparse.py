import GEOparse
import pandas as pd
import json
import os
import logging
import sys, os
module_dir = './'
sys.path.append(module_dir)
from src.constants import *
# It's assumed your 'helpers.py' file with 'get_geo_list' exists in the same directory.
from helpers import get_geo_list

# --- Configuration ---
# Set up basic logging to track progress and errors
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define paths for input and output
# GEO_DOWNLOAD_DIR = './downloads/geo_downloads/'
# METADATA_OUTPUT_DIR = './downloads/metadata/'
# COMBINED_DATA_OUTPUT_FILE = './downloads/combined_expression_data.csv'
# SOFT_PATH = './downloads/old_geo_downloads/'
# --- Helper Functions ---

def setup_directories():
    """Create output directories if they don't exist."""
    logging.info("Setting up output directories...")
    os.makedirs(GEO_DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(METADATA_OUTPUT_DIR, exist_ok=True)

def create_probe_to_gene_map(gpl):
    """
    Creates a mapping from probe IDs to gene symbols by searching for common
    gene symbol column names in the GPL table.
    """
    possible_gene_columns = ['ORF', 'Gene Symbol', 'GENE_SYMBOL', 'ILMN_Gene', 'Symbol']
    for col_name in possible_gene_columns:
        if col_name in gpl.table.columns:
            logging.info(f"Found gene symbol column: '{col_name}'")
            mapping_df = gpl.table[['ID', col_name]].dropna()
            return dict(zip(mapping_df['ID'], mapping_df[col_name]))
    logging.warning(f"Could not find a recognized gene symbol column in GPL {gpl.name}.")
    return None

def process_metadata(geo_accession, gse, gsm):
    """Filters and saves metadata for a given sample to the ./metadata folder."""
    excluded_keys = [
        'contact', 'date', '_id', 'taxid', 'data_', 'status', 'type', 
        'contributor', 'relation', 'geo_accession', 'row_count', 'organism', 
        'label', 'supplementary_file'
    ]
    filtered_study_meta = {k: v for k, v in gse.metadata.items() if not any(ex in k for ex in excluded_keys)}
    filtered_sample_meta = {k: v for k, v in gsm.metadata.items() if not any(ex in k for ex in excluded_keys)}
    final_metadata = {
        'study_id': geo_accession,
        'sample_id': gsm.name,
        'study_metadata': filtered_study_meta,
        'sample_metadata': filtered_sample_meta
    }
    output_path = os.path.join(METADATA_OUTPUT_DIR, f"{geo_accession}_{gsm.name}.json")
    try:
        with open(output_path, "w") as fp:
            json.dump(final_metadata, fp, indent=4)
    except Exception as e:
        logging.error(f"Failed to save metadata for {gsm.name}: {e}")

def process_sample_data(geo_accession, gsm, probe_to_gene_map):
    """
    Processes expression data for a sample and returns it as a DataFrame
    with a single, uniquely named column.
    """
    if gsm.table.empty:
        logging.warning(f"Sample {gsm.name} has no data table. Skipping.")
        return None

    df = gsm.table.copy()
    df['gene_symbol'] = df['ID_REF'].map(probe_to_gene_map)
    df.dropna(subset=['gene_symbol'], inplace=True)
    if df.empty:
        logging.warning(f"No genes could be mapped for sample {gsm.name}. Skipping.")
        return None

    # Filter for Arabidopsis thaliana genes (AtXgXXXXX format)
    at_gene_regex = r'^At[1-5MC]g\d{5}$'
    df = df[df['gene_symbol'].str.match(at_gene_regex, case=False)]
    if df.empty:
        logging.warning(f"No Arabidopsis thaliana genes found in {gsm.name}. Skipping.")
        return None

    df.set_index('gene_symbol', inplace=True)
    df = df[['VALUE']]
    
    # Handle duplicate genes by averaging their expression
    df = df.groupby(df.index).mean()
    
    # Create a unique name for the sample's data column
    unique_sample_id = f"{geo_accession}_{gsm.name}"
    df = df.rename(columns={'VALUE': unique_sample_id})
    
    return df

# --- Main Execution ---
# def main():
#     """
#     Main function to download, process, and combine GEO data into a single DataFrame
#     using a memory-efficient "collect and concatenate" pattern.
#     """
#     setup_directories()
    
#     try:
#         geo_list = get_geo_list('core_lists/data_addresses.csv')
#         master_gene_list_df = pd.read_csv('core_lists/genes_list.csv', index_col=0)
#     except FileNotFoundError as e:
#         logging.error(f"Core list file not found: {e}. Please ensure 'core_lists' directory is correct.")
#         return

#     # --- CHANGE 1: Prepare the master index and initialize an empty list ---
#     # We will collect all processed sample DataFrames in this list.
#     master_index_upper = master_gene_list_df.index.str.upper()
#     processed_samples_list = []
    
#     for geo_accession in geo_list:
#         logging.info(f"--- Processing study: {geo_accession} ---")
#         try:
#             gse = GEOparse.get_GEO(geo=geo_accession, destdir=GEO_DOWNLOAD_DIR, silent=True)
#         except Exception as e:
#             logging.error(f"Could not download or parse {geo_accession}. Error: {e}")
#             continue

#         try:
#             gpl_name = list(gse.gpls.keys())[0]
#             gpl = gse.gpls[gpl_name]
#             probe_to_gene_map = create_probe_to_gene_map(gpl)
#             if probe_to_gene_map is None:
#                 logging.error(f"Skipping study {geo_accession} due to missing gene map.")
#                 continue
#         except Exception as e:
#             logging.error(f"Error processing platform for {geo_accession}: {e}")
#             continue

#         for gsm_name, gsm in gse.gsms.items():
#             logging.info(f"Processing sample: {gsm_name}")
#             try:
#                 process_metadata(geo_accession, gse, gsm)
#                 processed_sample_df = process_sample_data(geo_accession, gsm, probe_to_gene_map)
                
#                 if processed_sample_df is not None and not processed_sample_df.empty:
#                     processed_sample_df.index = processed_sample_df.index.str.upper()
                    
#                     # --- CHANGE 2: Append the small DataFrame to the list ---
#                     # Instead of joining, we just add the result to our collection.
#                     # This is a very fast and low-memory operation.
#                     processed_samples_list.append(processed_sample_df)
#                     logging.info(f"Collected data for {gsm_name}.")

#             except Exception as e:
#                 logging.error(f"An unexpected error occurred while processing sample {gsm_name}: {e}")

#     # --- CHANGE 3: Assemble the final DataFrame after the loop ---
#     logging.info(f"--- All studies processed. Concatenating {len(processed_samples_list)} samples into the final DataFrame. ---")
    
#     if not processed_samples_list:
#         logging.warning("No samples were processed successfully. Output file will be empty.")
#         final_df = pd.DataFrame(index=master_index_upper)
#     else:
#         # This one-time concatenation is far more efficient than repeated joins.
#         all_samples_df = pd.concat(processed_samples_list, axis=1)
        
#         # Create an empty DataFrame with the master gene list index to ensure the final
#         # output conforms perfectly to that master list (same genes, same order).
#         final_df_scaffold = pd.DataFrame(index=master_index_upper)
        
#         # Join the concatenated data onto the master scaffold. This single join is fast.
#         final_df = final_df_scaffold.join(all_samples_df)

#     logging.info("Finalizing the combined DataFrame.")
#     final_df.dropna(axis=1, how='all', inplace=True)
#     final_df = final_df.loc[:, ~final_df.columns.duplicated()]
#     final_df.to_csv(COMBINED_DATA_OUTPUT_FILE)
#     logging.info(f"✅ Success! Combined DataFrame saved to '{COMBINED_DATA_OUTPUT_FILE}'")

def import_data():
    """
    Main function to download, process, and combine GEO data into a single DataFrame.
    """
    setup_directories()
    
    # Load the list of GEO studies and the master gene list
    try:
        geo_list = get_geo_list('core_lists/data_addresses.csv')
        master_gene_list_df = pd.read_csv('core_lists/genes_list.csv', index_col=0)
    except FileNotFoundError as e:
        logging.error(f"Core list file not found: {e}. Please ensure 'core_lists' directory is correct.")
        return

    # Initialize the final DataFrame with the index from the master gene list
    final_df = pd.DataFrame(index=master_gene_list_df.index)
    
    for geo_accession in geo_list:
        logging.info(f"--- Processing study: {geo_accession} ---")
        try:
            gse = GEOparse.get_GEO(geo=geo_accession, destdir=GEO_DOWNLOAD_DIR, silent=True)
        except Exception as e:
            logging.warning(f"Could not download or parse {geo_accession}. Error: {e}")
            try:
                gse = GEOparse.get_GEO(filepath=f"{SOFT_PATH}{geo_accession}_family.soft.gz")
            except Exception as ex:
                logging.error(f"Error whilt getting the local version of {geo_accession}. Error: {ex}")
                continue

        try:
            gpl_name = list(gse.gpls.keys())[0]
            gpl = gse.gpls[gpl_name]
            probe_to_gene_map = create_probe_to_gene_map(gpl)
            if probe_to_gene_map is None:
                logging.error(f"Skipping study {geo_accession} due to missing gene map, it is RNA seq.")
                continue
        except Exception as e:
            logging.error(f"Error processing platform for {geo_accession}: {e}")
            continue

        for gsm_name, gsm in gse.gsms.items():
            logging.info(f"Processing sample: {gsm_name}")
            try:
                # Goal 1: Store filtered metadata (unchanged)
                process_metadata(geo_accession, gse, gsm)
                
                # Goal 2: Process data and get a DataFrame to append
                processed_sample_df = process_sample_data(geo_accession, gsm, probe_to_gene_map)
                
                # If data was processed successfully, add it to the final DataFrame
                if processed_sample_df is not None and not processed_sample_df.empty:

                    processed_sample_df.index = processed_sample_df.index.str.upper()
                    processed_sample_df = processed_sample_df.loc[processed_sample_df.index.isin(final_df.index)]
                    final_df = pd.concat([final_df, processed_sample_df], axis=1)
                    # final_df = final_df.join(processed_sample_df)
                    logging.info(f"Appended data for {gsm_name} to the main DataFrame.")

            except Exception as e:
                logging.error(f"An unexpected error occurred while processing sample {gsm_name}: {e}")

    logging.info("--- All studies processed. Finalizing the combined DataFrame. ---")
    
    # After merging, some columns might be all NaN if they had no overlapping genes.
    # It's good practice to drop them.
    final_df.dropna(axis=1, how='all', inplace=True)

    # Save the final combined DataFrame to a single CSV file
    final_df.to_csv(COMBINED_DATA_OUTPUT_FILE)
    logging.info(f"✅ Success! Combined DataFrame saved to '{COMBINED_DATA_OUTPUT_FILE}'")

if __name__ == "__main__":
    import_data()

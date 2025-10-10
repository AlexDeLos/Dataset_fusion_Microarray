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
from data_importing.helpers import get_geo_list,mapping

# --- Configuration ---
# Set up basic logging to track progress and errors
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

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
    os.makedirs(PROCESSED_DATA_FOLDER, exist_ok=True)

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
    df['gene_symbol'] = df['gene_symbol'].map(lambda x: x.split("_")[0])
    if df.empty:
        logging.warning(f"No genes could be mapped for sample {gsm.name}. Skipping.")
        #TODO: a lot of the entries here already have proper mappings, so just let them through
        df = gsm.table.copy()
        df['gene_symbol'] = df['ID_REF'].map(lambda x: x.split("_")[0])
        # del df['ID_REF']
        # return None

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
def import_GEOparse():
    setup_directories()

    geo_list = get_geo_list('core_lists/data_addresses.csv')
    df = pd.read_csv('core_lists/genes_list.csv', index_col=0)
    df_index = pd.read_csv('core_lists/genes_list.csv', index_col=0)
    for number,geo in enumerate(geo_list):
        try:
            try:
                # gse =GEOparse.get_GEO(geo=geo, destdir=GEO_DOWNLOAD_DIR,silent=True)
                gse =GEOparse.get_GEO(filepath=f"{SOFT_PATH}{geo}_family.soft.gz")

            except FileNotFoundError as err:
                print(err)
            except EOFError as err:
                print(err)
            key = list(gse.gpls)
            key = key[0]
            gpl_table=gse.gpls[key].table
            try:
                probe_to_gene_map = dict(zip(gpl_table['ID'], gpl_table['ORF'].map(mapping)))
            except KeyError as error:
                print('---- ERROR MAKING MAP ----')
                print(error)
                continue
            for gsm_name, gsm in gse.gsms.items():
                in_df = gsm.table.copy()
                in_df['ID_REF'] = in_df['ID_REF'].map(probe_to_gene_map)
                # How to drop the numbers
                in_df = in_df.dropna()
                gsm_id = gsm_name+'_'+str(gse.metadata['type'])+'_'+str(number)
                in_df = in_df.rename(columns={'VALUE': gsm_id})
                in_df.set_index('ID_REF',inplace = True)
                # we take only the arabidopsis thaliana genges
                # in_df.drop(in_df[~in_df.index.str.match(r'^At[1-5MC]g\d{5}$', case=False)].index, inplace=True) # we take only the a
                in_df = in_df.loc[in_df.index.str.match(r'^At[1-5MC]g\d{5}$', case=False)]
                in_df = in_df.filter([gsm_id])

                if df.index.duplicated().any():
                    # print('Duplicates in df:', df[df.index.duplicated(keep=False)])

                    # Create a dictionary to keep track of the count of duplicates
                    duplicate_count = df.index.value_counts().to_dict()

                    # Group by the index column and calculate the mean of the values
                    df = df.groupby('ID_REF')[in_df.columns[0]].mean().reset_index() # assuming there 

                if in_df.index.duplicated().any():
                    # print('Duplicates in df:', in_df[in_df.index.duplicated(keep=False)])

                    # Create a dictionary to keep track of the count of duplicates
                    duplicate_count = in_df.index.value_counts().to_dict()

                    # Group by the index column and calculate the mean of the values
                    in_df = in_df.groupby('ID_REF')[in_df.columns[0]].mean().reset_index() # assuming there is only one value

                try:
                    if in_df.index.name != 'ID_REF':
                        in_df.set_index('ID_REF',inplace = True)
                    # Fill in all the nans
                    complete_in = pd.concat([df_index, in_df], axis=1)
                    # complete_in = complete_in.transform(lambda x: x.fillna(x.mean()))
                    df = pd.concat([df, complete_in], axis=1)
                    x = 0
                except Exception as error:
                    print(error)
                    pass
        


        except Exception as error:
            print(error)
            print('-----An error occured, probably an empty dataframe')
            
    df.to_csv(COMBINED_DATA_OUTPUT_FILE)
    # df = pd.DataFrame(index=df.index)
    logging.info(f"✅ Success! Combined DataFrame saved to '{COMBINED_DATA_OUTPUT_FILE}'")

def handle_duplicates(df_index,sample_df):

    if sample_df.index.duplicated().any():
        # print('Duplicates in df:', in_df[in_df.index.duplicated(keep=False)])

        # Create a dictionary to keep track of the count of duplicates
        duplicate_count = sample_df.index.value_counts().to_dict()

        # Group by the index column and calculate the mean of the values
        sample_df = sample_df.groupby('gene_symbol')[sample_df.columns[0]].mean().reset_index() # assuming there is only one value
    return sample_df


def import_data():
    """
    Main function to download, process, and combine GEO data into a single DataFrame.
    """
    setup_directories()
    
    # Load the list of GEO studies and the master gene list
    try:
        geo_list = Studies#get_geo_list('core_lists/data_addresses.csv')
        df_index = pd.read_csv('core_lists/genes_list.csv', index_col=0)
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

        gpl_names:list = list(gse.gpls.keys())# this run this for both gpls[0]
        for gpl_name in gpl_names:
            try:
                gpl_table = gse.gpls[gpl_name].table
                # probe_to_gene_map = create_probe_to_gene_map(gpl)
                probe_to_gene_map = dict(zip(gpl_table['ID'], gpl_table['ORF'].map(mapping)))
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
                    #TODO: check if GSM680342 is in the final df
                    processed_sample_df = process_sample_data(geo_accession, gsm, probe_to_gene_map)
                    
                    # If data was processed successfully, add it to the final DataFrame
                    if processed_sample_df is not None and not processed_sample_df.empty:
                        processed_sample_df.index = processed_sample_df.index.str.upper()
                        processed_sample_df = handle_duplicates(df_index,processed_sample_df)
                        complete_in = pd.concat([df_index, processed_sample_df], axis=1)
                        # complete_in = complete_in.transform(lambda x: x.fillna(x.mean()))
                        final_df = pd.concat([final_df, complete_in], axis=1)
                        # processed_sample_df = processed_sample_df.loc[processed_sample_df.index.isin(final_df.index)]
                        # final_df = pd.concat([final_df, processed_sample_df], axis=1)
                        # final_df = final_df.join(processed_sample_df)
                        logging.info(f"Appended data for {gsm_name} to the main DataFrame.")

                except Exception as e:
                    logging.error(f"An unexpected error occurred while processing sample {gsm_name}: {e}")

    logging.info("--- All studies processed. Finalizing the combined DataFrame. ---")
    
    # After merging, some columns might be all NaN if they had no overlapping genes.
    # It's good practice to drop them.
    # final_df.dropna(axis=1, how='all', inplace=True)

    # Save the final combined DataFrame to a single CSV file
    final_df.to_csv(COMBINED_DATA_OUTPUT_FILE)
    logging.info(f"✅ Success! Combined DataFrame saved to '{COMBINED_DATA_OUTPUT_FILE}'")

if __name__ == "__main__":
    # import_GEOparse()
    import_data()
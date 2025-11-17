from data_analisys.label_cluster_exploration import run_label_cluster_exploration
from data_analisys.diff_and_GSEA_pipeline import run_diff_exp_and_enrichment
from data_importing.import_GEOparse import import_data
from data_importing.data_norm_and_analisys import run_preprocessing
from meta_data_processing.label_generation import condense_labels
# RUN DATA IMPORTING
import_data()

#process and filter the data
run_preprocessing()
#RUN METADATA AND LABELING
condense_labels()

# RUN UMAP AND CLUSTER ANALISYS
run_label_cluster_exploration(0)
run_label_cluster_exploration(10)
run_label_cluster_exploration(15)

# RUN DIFF EXP and ENRICHMENT ANALISYS
run_diff_exp_and_enrichment()
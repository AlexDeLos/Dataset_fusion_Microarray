from src.data_analisys.label_cluster_exploration import run_label_cluster_exploration
from src.data_analisys.diff_and_GSEA_pipeline import run_diff_exp_and_enrichment
from src.data_importing.import_GEOparse import import_data
# RUN DATA IMPORTING
import_data()

#process and filter the data

# RUN UMAP AND CLUSTER ANALISYS
run_label_cluster_exploration()

# RUN DIFF EXP and ENRICHMENT ANALISYS
run_diff_exp_and_enrichment()
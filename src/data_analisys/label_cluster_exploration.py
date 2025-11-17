import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import silhouette_score
import numpy as np
import seaborn
import os
from itertools import compress
from sklearn.decomposition import PCA  # to apply PCA
import umap
import sys
import gc # Import the garbage collection module

module_dir = './'
sys.path.append(module_dir)
from src.data_analisys.utils.plot_utils import plot_tsne,plot_projection, plot_dendogram,plot_summary_scores_relative,plot_iteration_scores,plot_summary_scores_modified
from src.constants import *
from src.data_analisys.utils.cluster_exploration_utils import *
from src.data_analisys.utils.correlation import calculate_study_correlations

def apply_mask_to_maps(maps,mask):
    maps_ = maps.copy()
    for map in maps:
        map_ = list(compress(maps_[map], mask))
        maps_[map] = map_
    return maps_
def pca_variance(robust_df,path):
    pca = PCA()
    #
    # Determine transformed features
    #
    X_train_pca = pca.fit_transform(robust_df.T)
    #
    # Determine explained variance using explained_variance_ration_ attribute
    #
    exp_var_pca = pca.explained_variance_ratio_
    #
    # Cumulative sum of eigenvalues; This will be used to create step plot
    # for visualizing the variance explained by each principal component.
    #
    cum_sum_eigenvalues = np.cumsum(exp_var_pca)
    #
    # Create the visualization plot
    #
    plt.bar(range(0,len(exp_var_pca)), exp_var_pca, alpha=0.5, align='center', label='Individual explained variance')
    plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')
    plt.ylabel('Explained variance ratio')
    plt.xlabel('Principal component index')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(f'{path}/test.svg')
    plt.close()


def plot_var(df,path):
    plt.plot(np.sort(np.array(df.var())))
    # plt.plot(np.array(df.var()))
    plt.savefig(f'{path}/gene_var.svg')
    plt.close()
    
def run_label_cluster_exploration(fil=0):
    labels = load_labels_study(LABELS_PATH)
    calculate_study_correlations(labels)
    labels = keys_upper(labels)
    # Fuse the labels

    labels_types = ['TREATMENT','TISSUE','MEDIUM']
    labels_df = make_df_from_labels(labels, labels_types)
    del labels # Delete original labels dictionary as it's no longer needed


    labels_df['TREATMENT'] = labels_df['TREATMENT'].apply(lambda x: tuple(sorted(x)))
    labels_df['ID'] = labels_df.index
    labels_df.drop_duplicates(inplace=True)
    del labels_df['ID']
    labels_df['TREATMENT'] = labels_df['TREATMENT'].apply(lambda x: list(sorted(x)))
    outputr = f'{CLUSTER_EXPLORATION_FIGURES_DIR}{EXPERIMENT_NAME}'
    os.makedirs(outputr,exist_ok=True)
    labels_df.to_csv(f'{outputr}/labels.csv')

    full_scores = {}
    full_scores_sil = {}
    full_scores_w_study = {}
    TYPES = ['robust', 'standardized', '2_way_norm','study_corrected','imputed']#2_way_norm_og 'standardized+','2_way_norm_og',['2_way_norm','study_corrected','imputed']#
    # LINK_METHODS = ['single','complete','average','weighted','centroid','median']

    LINK_METHODS = ['complete']
    # for filt in [True,False]:
    for method in LINK_METHODS:
        experiment_name = f'{method}-{fil}'
        # experiment_name = f'norm_comp_15_BS_{method}'

        for type_ in TYPES:
            figure_out_path:str = f'{CLUSTER_EXPLORATION_FIGURES_DIR}{experiment_name}/{type_}'

            os.makedirs(figure_out_path,exist_ok=True)
            try:
                #'/home/alex/Documents/GitHub/Data_collection/df_final' they are the same as _og
                data_df = pd.read_csv(f'{PROCESSED_DATA_FOLDER}/{type_}.csv', index_col=0)
                if type_ =='2_way_norm_og':
                    data_df.columns = list(map(lambda x: f'{x.split('_')[-1]}_{x.split('_')[0]}',data_df.columns))
                # handel duplicate samples
                data_df = fuse_columns_by_sample(data_df)
            except FileNotFoundError:
                print(f"Warning: Data file for '{type_}' not found. Skipping.")
                continue
            except Exception as e:
                print(f"Error loading data for '{type_}': {e}")
                continue


            ## ----------------------------------------------------------------
            ## Filter out studies with fewer than 5 samples
            ## ----------------------------------------------------------------
            print(f"Processing data type: {type_}")
            print(f"Original shape: {data_df.shape}")
            # Get the study ID for each sample (column) and wrap in a pandas Series
            studies_series = pd.Series(get_studies(data_df), index=data_df.columns)

            # Count the number of samples per study
            study_counts = studies_series.value_counts()

            # Identify studies that have 5 or more samples
            studies_to_keep = study_counts[study_counts >= fil].index

            # If no studies meet the criteria, skip this entire dataframe type
            if studies_to_keep.empty:
                print(f"--> Skipping '{type_}' as no studies have 5 or more samples.")
                del data_df, studies_series, study_counts
                gc.collect()
                continue # Move to the next `type_`

            # Create a boolean mask and apply it to keep only samples from the desired studies
            mask = studies_series.isin(studies_to_keep)
            data_df = data_df.loc[:, mask]

            print(f"Filtered shape (studies with >= {fil} samples): {data_df.shape}")
            ## ----------------------------------------------------------------
            ## End of new filtering logic
            ## ----------------------------------------------------------------


            samples = get_samples(data_df)
            studies = get_studies(data_df)

            # Save the labels of only the samples I use
            labels_df[labels_df.index.isin(samples)]
            labels_df.to_csv(f'{figure_out_path}/labels.csv')
            #TODO: investigate the labels that have no sample linked to it

            maps = get_label_map_new(data_df,labels_df)

            linkage_method = method
            number_of_clusters = 15

            pca = PCA(n_components = 50)
            pca_variance(data_df,path=figure_out_path)
            plot_var(data_df.T,path=figure_out_path)
            pca.fit(data_df.T)
            data_pca = pca.transform(data_df.T)
            data_pca = pd.DataFrame(data_pca)
            maps = add_map(maps,hierarchy.fcluster(hierarchy.linkage(data_pca, method=linkage_method),t=number_of_clusters,criterion='maxclust'),'PCA')
            plot_dendogram(data_pca,linkage_method,number_of_clusters,figure_out_path,name='50D PCA')
            del data_pca # Free memory from PCA results


            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5)
            embedding = reducer.fit_transform(data_df.T.to_numpy())
            plot_dendogram(embedding,linkage_method,number_of_clusters,figure_out_path,name='Umap 2D embeding')
            maps = add_map(maps,hierarchy.fcluster(hierarchy.linkage(embedding, method=linkage_method),t=number_of_clusters,criterion='maxclust'),'emb_2D')

            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5,n_components=50)
            embedding_50 = reducer.fit_transform(data_df.T.to_numpy())
            new_temp = hierarchy.linkage(embedding_50, method=linkage_method)

            plot_dendogram(embedding_50,linkage_method,number_of_clusters,figure_out_path,name='Umap 50D embeding')
            maps = add_map(maps,hierarchy.fcluster(new_temp,t=number_of_clusters,criterion='maxclust'),'emb_50D')
            del embedding_50 # Free memory from 50D UMAP embedding
            del new_temp # Free memory from linkage matrix

            temp_old = hierarchy.linkage(data_df.T.to_numpy(), method=linkage_method)
            cluster = hierarchy.fcluster(temp_old,t=number_of_clusters,criterion='maxclust') #len(np.unique(study_map))
            plot_dendogram(temp_old,linkage_method,number_of_clusters,figure_out_path,name='No dim reduction')
            maps = add_map(maps,cluster,'clusters')

            def get_optimal_clusters(data):
                linkage = hierarchy.linkage(data.T.to_numpy(), method=linkage_method)
                cluster = hierarchy.fcluster(linkage,t=0.1,criterion='inconsistent')
                # calculate the fitnes of that number
                del linkage # Free memory from linkage matrix inside function
                return cluster
            cluster = get_optimal_clusters(data_df)
            maps = add_map(maps,cluster,'clusters_use')
            for sample in maps:

            # # maps['clusters_1'] = cluster_1
                # del maps['study']
                del maps[sample]['clusters']
                del maps[sample]['clusters_use']
                del maps[sample]['MEDIUM']
                # del maps['tissue']
                del maps[sample]['PCA']
                del maps[sample]['emb_2D']
                del maps[sample]['emb_50D']
                # del maps['agar']

            # Plot pie charts
            plot_pie_chart(maps,figure_out_path)

            #! temp clutering on low dim:

            #TODO: scores for treatment should be calculated by tissue
            # Let's start with the leaf tissue, as it is the biggest
            exmaple_sample = list(maps.keys())[0]
            scores = list(map(lambda x: {x:adjusted_rand_score(
                hierarchy.fcluster(temp_old,t=len(np.unique(get_map(maps,x))),criterion='maxclust'),get_map(maps,x)
                )}, maps[exmaple_sample]))

            del temp_old # Free memory from the large linkage matrix

            scores_with_study = list(map(lambda x: {x:adjusted_rand_score(get_map(maps,'study'),get_map(maps,x))}, maps[exmaple_sample]))

            sil_scores = list(map(lambda x: {x:silhouette_score(data_df.T.to_numpy(),labels=get_map(maps,x))}, maps[exmaple_sample]))

            #TODO: re-add this
            tissue_scores = []
            # for tissue in ['leaf']:
            #     mask_tissue = np.array([x ==tissue for x in get_map(maps,'TISSUE')])
            #     tiss_df = data_df.loc[:, mask_tissue]
            #     temp_treatment = hierarchy.linkage(tiss_df.T.to_numpy(), method=linkage_method)
            #     maps_tissue = apply_mask_to_maps(maps,mask_tissue)#TODO: this needs to be fixed to the new maps object

            #     tissue_scores.append(adjusted_rand_score(hierarchy.fcluster(temp_treatment,t=len(np.unique(maps_tissue['TREATMENT'])),criterion='maxclust'),maps_tissue['TREATMENT']))

            for i,key in enumerate(maps[exmaple_sample]):
                plot_tsne(
                    df=data_df.T,
                    markers=get_map(maps,key),
                    colors=get_map(maps,key),
                    save_path= f'{figure_out_path}/tsne',
                    title=None,#f'TSNE projection for {key}',
                    name=key,
                    legend=True)
                plot_projection(embedding,
                    markers=get_map(maps,key),
                    colors=get_map(maps,key),
                    title=None,#f'TSNE projection for {key}',
                    name=key,
                    legend=True,  # Shows legend for color groups
                    save_path=f'{figure_out_path}/Umaps')

            del embedding # Free memory from 2D UMAP embedding

            def get_color_df(maps,scores)->pd.DataFrame:
                dic = {}
                col = 0
                for key in maps[exmaple_sample]:
                    int_map = to_int(get_map(maps,key),name=key,path=figure_out_path)
                    palette = seaborn.color_palette("husl", len(int_map)+1)  # Choose a palette
                    col_colors = [palette[i] for i in int_map]  # Map cluster IDs to colors
                    dic[f'{key,"{:0.3e}".format(scores[col][key])}']= col_colors
                    col = col +1
                return pd.DataFrame(dic, index=data_df.columns)
            color_df = get_color_df(maps,scores)

            # plot_heat_map(data_df,figure_out_path,cluster=False,col_cluster=False,typ='png',title = 'overview overview study',log_norm=True, col=color_df, name='general_overview_study')

            print(f'rand indx {type_}: {scores}')
            print(f'sil {type_}: {sil_scores}')

            for i,el in enumerate(scores):
                key = list(el.keys())[0]
                if 'cluster' in key:
                    pass
                else:
                    full_scores[f'{type_} val {key}'] = scores[i][key]
            for i,el in enumerate(sil_scores):
                key = list(el.keys())[0]
                if 'cluster' in key:
                    pass
                else:
                    full_scores_sil[f'{type_} val {key}'] = sil_scores[i][key]
            for i,el in enumerate(scores_with_study):
                key = list(el.keys())[0]
                full_scores_w_study[f'{type_} val {key}'] = scores_with_study[i][key]

            ## ----------------------------------------------------------------
            ## START: REFACTORED PLOTTING LOGIC
            ## ----------------------------------------------------------------
            # Consolidate scores into single dictionaries
            current_sil_scores = {key: value for score_dict in sil_scores for key, value in score_dict.items()}
            current_rand_scores = {key: value for score_dict in scores for key, value in score_dict.items()}

            # Call the new plotting function for Silhouette scores
            plot_iteration_scores(
                scores_dict=current_sil_scores,
                y_label='Silhouette Score',
                title=f'Silhouette Scores for {type_} ({method} linkage)',
                file_name='silhouette_scores.svg',
                output_path=figure_out_path
            )
            # plot_iteration_scores_modified(
            #     scores_dict=current_sil_scores,
            #     y_label='Silhouette Score Comparison',
            #     title=f'Silhouette Scores for {type_} ({method} linkage)',
            #     file_name='silhouette_scores_compare.svg',
            #     output_path=figure_out_path
            # )

            # Call the new plotting function for Rand Index scores
            plot_iteration_scores(
                scores_dict=current_rand_scores,
                y_label='Rand Index Score',
                title=f'Rand Index Scores for {type_} ({method} linkage)',
                file_name='rand_index_scores.svg',
                output_path=figure_out_path
            )
            # plot_iteration_scores_modified(
            #     scores_dict=current_rand_scores,
            #     y_label='Rand Index Score Comparison',
            #     title=f'Rand Index Scores for {type_} ({method} linkage)',
            #     file_name='rand_index_scores_compare.svg',
            #     output_path=figure_out_path
            # )
            ## ----------------------------------------------------------------
            ## END: REFACTORED PLOTTING LOGIC
            ## ----------------------------------------------------------------

            # Clean up all large objects at the end of the loop iteration
            del data_df, maps, scores, sil_scores, scores_with_study, color_df, cluster
            gc.collect() # Trigger garbage collection to free up memory

    # --- FINAL SUMMARY PLOTTING CALLS ---
    output_dir = f'{CLUSTER_EXPLORATION_FIGURES_DIR}{EXPERIMENT_NAME}/{fil}/'

    # Call the summary plotting function for each score type
    plot_summary_scores_modified(
        scores_dict=full_scores,
        title='Rand Index Score',
        file_name='barRandInd.svg',
        output_dir=output_dir
    )
    plot_summary_scores_relative(
        scores_dict=full_scores,
        title='Rand Index Score Ralative to Study Silhouette Rand Index Score',
        file_name='compBarRandInd.svg',
        output_dir=output_dir
    )

    plot_summary_scores_modified(
        scores_dict=full_scores_sil,
        title='Silhouette Score',
        file_name='barSilhouette.svg',
        output_dir=output_dir
    )
    plot_summary_scores_relative(
        scores_dict=full_scores_sil,
        title='Silhouette Score Ralative to Study Silhouette Score',
        file_name='compBarSilhouette.svg',
        output_dir=output_dir
    )

    print("DONE WITH CLUSTER EXPLORATION")
if __name__ == '__main__':
    # run_label_cluster_exploration()
    run_label_cluster_exploration(15)
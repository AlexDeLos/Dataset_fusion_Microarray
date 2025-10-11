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
module_dir = './'
sys.path.append(module_dir)
from src.data_analisys.utils.plot_utils import plot_tsne,get_Umap_3,plot_heat_map, plot_dendogram
from src.constants import *
from src.data_analisys.utils.cluster_exploration_utils import *

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

# def get_study(sample: str):
#     return sample.split('_')[0]

def run_label_cluster_exploration():
    labels = load_labels_study(LABELS_PATH)
    labels = keys_upper(labels) 
    # Fuse the labels

    labels_types = ['TREATMENT','TISSUE']#,'MEDIUM']
    labels_df = make_df_from_labels(labels, labels_types)
    full_scores = {}
    full_scores_sil = {}
    full_scores_w_study = {}
    TYPES = ['robust', 'standardized', 'robust+', 'standardized+','2_way_norm','study_corrected','imputed']
    # LINK_METHODS = ['single','complete','average','weighted','centroid','median']

    LINK_METHODS = ['complete']
    for method in LINK_METHODS:
        experiment_name = f'{EXPERIMENT_NAME}-{method}'
        # experiment_name = f'norm_comp_15_BS_{method}'

        for type_ in TYPES:
            figure_out_path:str = f'{CLUSTER_EXPLORATION_FIGURES_DIR}{experiment_name}-{type_}'

            os.makedirs(figure_out_path,exist_ok=True)
            try:
                data_df = pd.read_csv(f'{PROCESSED_DATA_FOLDER}/{type_}.csv', index_col=0)
                # handel duplicate samples
                data_df = fuse_columns_by_sample(data_df)
            except:
                raise ValueError

            # low_sample_study_filter = list(map(lambda x: int(x) if x > 10 else -1,np.unique(np.array(data_df.T.index.map(lambda x: int(x.split('_')[-1]))),return_counts=True)[1]))
            # TODO filter out the smaller dataset?
            # Apply filter to df
            # robust_df.T = robust_df.T.loc[:,robust_df.T.index.map(lambda x: int(x.split('_')[-1]))wait]

            samples = get_samples(data_df)
            studies = get_studies(data_df)

            # Save the labels of only the samples I use
            labels_df[labels_df.index.isin(samples)]
            labels_df['TREATMENT'] = labels_df['TREATMENT'].apply(lambda x: tuple(sorted(x)))
            labels_df['ID'] = labels_df.index
            labels_df.drop_duplicates(inplace=True)
            del labels_df['ID']
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


            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5)
            embedding = reducer.fit_transform(data_df.T.to_numpy())
            plot_dendogram(embedding,linkage_method,number_of_clusters,figure_out_path,name='Umap 2D embeding')

            maps = add_map(maps,hierarchy.fcluster(hierarchy.linkage(embedding, method=linkage_method),t=number_of_clusters,criterion='maxclust'),'emb_2D')


            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5,n_components=50)
            embedding_50 = reducer.fit_transform(data_df.T.to_numpy())
            new_temp = hierarchy.linkage(embedding_50, method=linkage_method)

            plot_dendogram(embedding_50,linkage_method,number_of_clusters,figure_out_path,name='Umap 50D embeding')

            maps = add_map(maps,hierarchy.fcluster(new_temp,t=number_of_clusters,criterion='maxclust'),'emb_50D')


            temp_old = hierarchy.linkage(data_df.T.to_numpy(), method=linkage_method)
            cluster = hierarchy.fcluster(temp_old,t=number_of_clusters,criterion='maxclust') #len(np.unique(study_map))
            plot_dendogram(temp_old,linkage_method,number_of_clusters,figure_out_path,name='No dim reduction')
            maps = add_map(maps,cluster,'clusters')

            def get_optimal_clusters(data):
                linkage = hierarchy.linkage(data.T.to_numpy(), method=linkage_method)
                cluster = hierarchy.fcluster(linkage,t=0.1,criterion='inconsistent')
                # calculate the fitnes of that number
                return cluster
            cluster = get_optimal_clusters(data_df)
            maps = add_map(maps,cluster,'clusters_use')
            for sample in maps:

            # # maps['clusters_1'] = cluster_1
                # del maps['study']
                del maps[sample]['clusters']
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
            scores_with_study = list(map(lambda x: {x:adjusted_rand_score(get_map(maps,'study'),get_map(maps,x))}, maps[exmaple_sample]))
            
            sil_scores = list(map(lambda x: {x:silhouette_score(data_df.T.to_numpy(),labels=get_map(maps,x))}, maps[exmaple_sample]))

            #TODO: re-add this
            # for tissue in ['leaf']:
            #     mask_tissue = np.array([x ==tissue for x in get_map(maps,'TISSUE')])
            #     tiss_df = data_df.loc[:, mask_tissue]
            #     temp_treatment = hierarchy.linkage(tiss_df.T.to_numpy(), method=linkage_method)
            #     maps_tissue = apply_mask_to_maps(maps,mask_tissue)#TODO: this needs to be fixed to the new maps object

            #     scores.append({f'treatment_{tissue}': adjusted_rand_score(hierarchy.fcluster(temp_treatment,t=len(np.unique(maps_tissue['TREATMENT'])),criterion='maxclust'),maps_tissue['TREATMENT'])})
                
            #     sil_scores.append({f'treatment_{tissue}':silhouette_score(tiss_df.T.to_numpy(),labels=maps_tissue['TREATMENT'])})
            
            # Perform UMAP transformation
            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5)
            embedding = reducer.fit_transform(data_df.T.to_numpy())
            for i,key in enumerate(maps[exmaple_sample]):
                plot_tsne(
                    df=data_df.T,
                    markers=get_map(maps,key),
                    colors=get_map(maps,key),
                    save_path= f'{figure_out_path}/tsne',
                    title=f'TSNE of {key} using {method} linkage',
                    name=key,
                    legend=True,  # Shows legend for color groups
                )
                get_Umap_3(embedding,name=key,colour_map=get_map(maps,key),marker_map=None, title =f'Umap of {key} using {method} linkage', save_loc=f'{figure_out_path}/Umaps')


            def get_color_df(maps,scores):
                dic = {}
                col = 0
                for key in maps[exmaple_sample]:
                    int_map = to_int(get_map(maps,key),name=key,path=figure_out_path)
                    palette = seaborn.color_palette("husl", len(int_map)+1)  # Choose a palette
                    col_colors = [palette[i] for i in int_map]  # Map cluster IDs to colors
                    dic[f'{key,'{:0.3e}'.format(scores[col][key])}']= col_colors
                    col = col +1
                return pd.DataFrame(dic, index=data_df.columns)
            color_df = get_color_df(maps,scores)
            
            plot_heat_map(data_df,figure_out_path,cluster=False,col_cluster=False,typ='png',title = 'overview overview study',log_norm=True, col=color_df, name='general_overview_study')
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
            # for i,k in enumerate(maps):
            #     if 'cluster' in k:
        # study_map = list(map(get_study,df_impute.columns))
        # # raise ValueError("Running in the cluster")
        # d = dict([(y,x+1) for x,y in enumerate(sorted(set(study_map)))])
        # batches = []
        # for el in study_map:
        #     batches.append(d[el])
        #     #         pass
        #     #     else:
        #     #         full_scores[f'{type_} val {k}'] = scores[i][k]
        #     # for i,k in enumerate(maps):
        #     #     if 'cluster' in k:
        #     #         pass
        #     #     else:
        #     #         full_scores_sil[f'{type_} val {k}'] = sil_scores[i][k]
        
    plt.bar(range(len(full_scores)), list(full_scores.values()), align='center')
    plt.xticks(range(len(full_scores)), list(full_scores.keys()),rotation = 90)
    plt.tight_layout()
    plt.title('rand ind score')
    plt.savefig(f'{CLUSTER_EXPLORATION_FIGURES_DIR}/{EXPERIMENT_NAME}/bar_rand_ind.svg')
    plt.close()

    plt.bar(range(len(full_scores_sil)), list(full_scores_sil.values()), align='center')
    plt.xticks(range(len(full_scores_sil)), list(full_scores_sil.keys()),rotation = 90)
    plt.tight_layout()
    plt.title('silhouette score')
    plt.savefig(f'{CLUSTER_EXPLORATION_FIGURES_DIR}/{EXPERIMENT_NAME}/bar_silhouette.svg')
    plt.close()


    plt.bar(range(len(full_scores_w_study)), list(full_scores_w_study.values()), align='center')
    plt.xticks(range(len(full_scores_w_study)), list(full_scores_w_study.keys()),rotation = 90)
    plt.tight_layout()
    plt.title('rand ind score study')
    plt.savefig(f'{CLUSTER_EXPLORATION_FIGURES_DIR}/{EXPERIMENT_NAME}/bar_rand_ind_study.svg')
    plt.close()

if __name__ =='__main__':
    run_label_cluster_exploration()
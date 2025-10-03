import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from utils.plot_utils import plot_tsne,get_Umap_3,plot_heat_map, plot_dendogram
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import silhouette_score
import numpy as np
import seaborn
import os
from itertools import compress
from utils.cluster_exploration_utils import *
from sklearn.decomposition import PCA  # to apply PCA
import umap

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

save_dir:str = '/home/alex/Documents/GitHub/Data_collection/df_final'


def run_label_cluster_exploration(save_dir:str= '/home/alex/Documents/GitHub/Data_collection/df_final', labels_path:str='/home/alex/Documents/GitHub/meta_data/in_use_labels/venice'):
    labels = load_labels_study(labels_path)
    labels = keys_upper(labels) 
    # Fuse the labels

    labels_types = ['TREATMENT','TISSUE','MEDIUM']
    df = make_df_from_labels(labels, labels_types)

    SAVE:bool = True
    TOP_N:int = -1
    full_scores = {}
    full_scores_sil = {}
    full_scores_w_study = {}
    TYPES = ['2_way_norm','corrected']
    # LINK_METHODS = ['single','complete','average','weighted','centroid','median']

    LINK_METHODS = ['complete']
    for method in LINK_METHODS:
        experiment_name = f'1.0-{method}'
        # experiment_name = f'norm_comp_15_BS_{method}'
        for type_ in TYPES:
            figure_out_path:str = f'./outputs/data_outputs/{experiment_name}-{type_}'

            os.makedirs(figure_out_path,exist_ok=True)
            try:
                robust_df = pd.read_csv(f'{save_dir}/{type_}.csv', index_col=0)
            except:
                # raise ValueError
                data = pd.read_csv(f'{save_dir}/corrected.csv', index_col=0)
                mat = data.to_numpy()
                q75, q25 = np.percentile(mat, [75 ,25],axis=1,keepdims=True)
                iqr = q75 - q25
                norm = (mat - np.median(mat,axis=1,keepdims=True))/iqr
                robust_df = pd.DataFrame(norm, columns=data.columns,index=data.index)

                robust_df.to_csv(f'{save_dir}/2_way_nrom.csv')

            #!single study => 155 is GSE41935
            #!single study => 116 is GSE34188
            #!single study => 40 is GSE40061
            #!single study => 64 is GSE27550


            #top 5 biggest-> 155,116,7,196,65
            # 196 -> GSE201609
            # 7 ->GSE46205
            # 65 -> GSE27548
            # 155 is GSE41935
            # 116 is GSE34188

            # 64 is GSE27550
            # 3 is GSE76827

            # 5 is GSE5620 	
            # 18 is GSE110079
            # 172 is GSE60960

            # 2 is GSE5624
            # 4 is GSE5622
            # 93 is GSE126373

            # TODO: check labels for these studies
            # 117 is GSE5628
            # 145 is GSE63128
            # robust_df = robust_df.loc[:, robust_df.columns.str.contains('_(116|196|7|65|155|64|5|18|172|2|3|4|93|117|145)$')]
            # robust_df = robust_df.loc[:, robust_df.columns.str.contains('_(116|196|7|65|64|5|18|172|2|3|4|93|117|145|155)$')]#155 missing labels

            low_sample_study_filter = list(map(lambda x: int(x) if x > 10 else -1,np.unique(np.array(robust_df.T.index.map(lambda x: int(x.split('_')[-1]))),return_counts=True)[1]))
            # TODO filter out the smaller dataset?
            # Apply filter to df
            # robust_df.T = robust_df.T.loc[:,robust_df.T.index.map(lambda x: int(x.split('_')[-1]))wait]
            
            study_map = list(robust_df.T.index.map(lambda x: int(x.split('_')[-1])))

            sample_temp = list(robust_df.T.index.map(lambda x: x.split('_')[0]))
            study_sample_map = zip(study_map,sample_temp)

            new_col = list(map(lambda x : x.split('_')[0],robust_df.columns))
            robust_df.columns = new_col



            #DEBUG CODE TODO: Delete
            seen = set()
            # A list to store duplicates found in the input list
            duplicates = []

            # Iterate over each element in the list
            for i in list(df.index):
                if i in seen:
                    duplicates.append(i)
                else:
                    seen.add(i)


            # Save the labels of only the samples I use
            df[df.index.isin(new_col)]
            df.to_csv(f'{figure_out_path}/labels.csv')
            get_biggest_studies(study_map,new_col,15)
            #SORT DF
            listings = list(robust_df.columns)
            listings.sort()
            robust_df = robust_df[listings]
            sample_index = list(robust_df.columns)

            #Fix study map order
            sorted_study = sorted(study_sample_map,key=lambda k: k[1])
            study_map = list(map( lambda g: g[0],sorted_study))

            study_map,seen = to_int_(study_map)
            labels_array = []
            samples_array = []
            for key in labels:
                for sample in labels[key]:
                    samples_array.append(sample)
                    save_dic = labels[key][sample]
                    save_dic['sample_id'] = sample
                    labels_array.append(labels[key][sample])
            newlist = sorted(labels_array, key=lambda d: d['sample_id'])
            maps = get_label_map(labels=newlist, sample_index=sample_index,figure_out_path=figure_out_path,labels_types=labels_types)
            maps['study'] = study_map

            missing_labels = sum(x is None for x in maps[labels_types[0]])
            get_incomplete_studies(maps,sample_index,labels_types[0])
            # Clear the fact that there might be no labels
            mask = np.array([x is not None for x in maps[labels_types[0]]])

            maps = apply_mask_to_maps(maps,mask)
            robust_df = robust_df.loc[:, mask]


            linkage_method = method
            number_of_clusters = 15

            pca = PCA(n_components = 50)
            pca_variance(robust_df,path=figure_out_path)
            plot_var(robust_df.T,path=figure_out_path)
            pca.fit(robust_df.T)
            data_pca = pca.transform(robust_df.T)
            data_pca = pd.DataFrame(data_pca)
            maps['PCA'] = hierarchy.fcluster(hierarchy.linkage(data_pca, method=linkage_method),t=number_of_clusters,criterion='maxclust')
            plot_dendogram(data_pca,linkage_method,number_of_clusters,figure_out_path,name='50D PCA')


            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5)
            embedding = reducer.fit_transform(robust_df.T.to_numpy())
            plot_dendogram(embedding,linkage_method,number_of_clusters,figure_out_path,name='Umap 2D embeding')

            maps['emb_2D'] = hierarchy.fcluster(hierarchy.linkage(embedding, method=linkage_method),t=number_of_clusters,criterion='maxclust')


            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5,n_components=50)
            embedding_50 = reducer.fit_transform(robust_df.T.to_numpy())
            new_temp = hierarchy.linkage(embedding_50, method=linkage_method)

            plot_dendogram(embedding_50,linkage_method,number_of_clusters,figure_out_path,name='Umap 50D embeding')

            maps['emb_50D'] = hierarchy.fcluster(new_temp,t=number_of_clusters,criterion='maxclust')


            temp_old = hierarchy.linkage(robust_df.T.to_numpy(), method=linkage_method)
            cluster = hierarchy.fcluster(temp_old,t=number_of_clusters,criterion='maxclust') #len(np.unique(study_map))
            plot_dendogram(temp_old,linkage_method,number_of_clusters,figure_out_path,name='No dim reduction')
            maps['clusters'] = cluster

            def get_optimal_clusters(data):
                linkage = hierarchy.linkage(data.T.to_numpy(), method=linkage_method)
                cluster = hierarchy.fcluster(linkage,t=0.1,criterion='inconsistent')
                # calculate the fitnes of that number
                return cluster
            cluster = get_optimal_clusters(robust_df)
            maps['clusters_use'] = cluster


            for key in maps:
                maps[key] = list(map(lambda x: str(x), maps[key]))
            # maps['clusters_1'] = cluster_1
            # del maps['study']
            del maps['clusters']
            # del maps['tissue']
            del maps['PCA']
            del maps['emb_2D']
            del maps['emb_50D']
            # del maps['agar']

            # Plot pie charts
            plot_pie_chart(maps,figure_out_path)

            #! temp clutering on low dim:

            #TODO: scores for treatment should be calculated by tissue
            # Let's start with the leaf tissue, as it is the biggest

            scores = list(map(lambda x: {x:adjusted_rand_score(
                hierarchy.fcluster(temp_old,t=len(np.unique(maps[x])),criterion='maxclust'),maps[x]
                )}, maps))
            scores_with_study = list(map(lambda x: {x:adjusted_rand_score(maps['study'],maps[x])}, maps))
            
            sil_scores = list(map(lambda x: {x:silhouette_score(robust_df.T.to_numpy(),labels=maps[x])}, maps))

            for tissue in ['leaf']:
                mask_tissue = np.array([x ==tissue for x in maps['TISSUE']])
                tiss_df = robust_df.loc[:, mask_tissue]
                temp_treatment = hierarchy.linkage(tiss_df.T.to_numpy(), method=linkage_method)
                maps_tissue = apply_mask_to_maps(maps,mask_tissue)

                scores.append({f'treatment_{tissue}': adjusted_rand_score(hierarchy.fcluster(temp_treatment,t=len(np.unique(maps_tissue['TREATMENT'])),criterion='maxclust'),maps_tissue['TREATMENT'])})
                
                sil_scores.append({f'treatment_{tissue}':silhouette_score(tiss_df.T.to_numpy(),labels=maps_tissue['TREATMENT'])})
            
            # Perform UMAP transformation
            reducer = umap.UMAP(n_epochs=200,n_neighbors=500, min_dist= 0.5)
            embedding = reducer.fit_transform(robust_df.T.to_numpy())
            for i,key in enumerate(maps):
                plot_tsne(
                    df=robust_df.T,
                    markers=maps[key],#ierarchy.fcluster(temp_old,t=len(np.unique(maps[key])),criterion='maxclust'),
                    colors=maps[key],
                    save_path= f'{figure_out_path}/tsne',
                    title=f'TSNE of {key} using {method} linkage',
                    name=key,
                    legend=True,  # Shows legend for color groups
                )
                get_Umap_3(embedding,name=key,colour_map=maps[key],marker_map=None, title =f'Umap of {key} using {method} linkage', save_loc=f'{figure_out_path}/Umaps')


            def get_color_df(maps,scores):
                dic = {}
                col = 0
                for m in maps:
                    int_map = to_int(maps[m],name=m,path=figure_out_path)
                    palette = seaborn.color_palette("husl", len(int_map)+1)  # Choose a palette
                    col_colors = [palette[i] for i in int_map]  # Map cluster IDs to colors
                    dic[f'{m,'{:0.3e}'.format(scores[col][m])}']= col_colors
                    col = col +1
                return pd.DataFrame(dic, index=robust_df.columns)
            color_df = get_color_df(maps,scores)
            
            plot_heat_map(robust_df,figure_out_path,cluster=False,col_cluster=False,typ='png',title = 'overview overview study',log_norm=True, col=color_df, name='general_overview_study')
            test = plot_heat_map(robust_df,figure_out_path,cluster=False,col_cluster=True,typ='png',title = 'overview overview study',log_norm=True, col=color_df, name='general_overview')
            # del scores['key']
            # del scores['key']
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
            #         pass
            #     else:
            #         full_scores[f'{type_} val {k}'] = scores[i][k]
            # for i,k in enumerate(maps):
            #     if 'cluster' in k:
            #         pass
            #     else:
            #         full_scores_sil[f'{type_} val {k}'] = sil_scores[i][k]
        
    plt.bar(range(len(full_scores)), list(full_scores.values()), align='center')
    plt.xticks(range(len(full_scores)), list(full_scores.keys()),rotation = 90)
    plt.tight_layout()
    plt.title('rand ind score')
    plt.savefig(f'{figure_out_path}/bar_rand_ind.svg')
    plt.close()

    plt.bar(range(len(full_scores_sil)), list(full_scores_sil.values()), align='center')
    plt.xticks(range(len(full_scores_sil)), list(full_scores_sil.keys()),rotation = 90)
    plt.tight_layout()
    plt.title('silhouette score')
    plt.savefig(f'{figure_out_path}/bar_silhouette.svg')
    plt.close()


    plt.bar(range(len(full_scores_w_study)), list(full_scores_w_study.values()), align='center')
    plt.xticks(range(len(full_scores_w_study)), list(full_scores_w_study.keys()),rotation = 90)
    plt.tight_layout()
    plt.title('rand ind score study')
    plt.savefig(f'{figure_out_path}/bar_rand_ind_study.svg')
    plt.close()

if __name__ =='__main__':
    run_label_cluster_exploration()
import json
import pandas as pd
import numpy as np
import os
import math
import matplotlib.pyplot as plt
import re

def to_int(x:list,name:str,path:str):
    new_x = []
    seen = {}
    floor = 1
    for i in x:
        if type(i) == list:
            # we need to find a way to deal with lists
            i.sort()
            i = str(i)
        if type(i) == dict:
            i = '?'
        if i in seen:
            new_x.append(seen[i])
        elif i is None:
            new_x.append(0)
        else:
            new_x.append(floor)
            seen[i]= floor
            floor = floor+1
    save_json(object=seen,name=f'{name}_dict',path=path)
    return new_x

def to_int_(x:list):
    new_x = []
    seen = {}
    floor = 1
    for i in x:
        if type(i) == list:
            # we need to find a way to deal with lists
            i.sort()
            i = str(i)
        if type(i) == dict:
            i = '?'
        if i in seen:
            new_x.append(seen[i])
        elif i is None:
            new_x.append(0)
        else:
            new_x.append(floor)
            seen[i]= floor
            floor = floor+1
    return new_x,seen

def get_label_map(labels,sample_index,figure_out_path,labels_types):
    labels_map = {}
    for key in labels_types:
        control_map = [None] * len(sample_index)
        for sample in labels:
            try:
                control_map[sample_index.index(sample['sample_id'].upper())] = sample[key] #[Control]
            except ValueError as e:
                try:
                    control_map[sample_index.index(sample['sample_id'].upper())] = '?'
                except:
                    pass
        labels_map[key] = control_map#to_int(control_map,name=key,path=figure_out_path)
    return labels_map


def flatten_extend(matrix):
    flat_list = []
    for row in matrix:
        if type(row) is list:
            flat_list.extend(row)
        else:
            flat_list = matrix
    return flat_list

def sanatize_labels(labels):
    for study in labels:
        for sample in labels[study]:
            if labels[study][sample]['treatment'] is None:
                labels[study][sample]['treatment'] = [None]
            
            labels[study][sample]['treatment'] = flatten_extend(labels[study][sample]['treatment'])
            for index in range(len(labels[study][sample]['treatment'])):
                st = False
                if type(labels[study][sample]['treatment']) is str:
                    st= True
                    labels[study][sample]['treatment'] = [labels[study][sample]['treatment']]
                treatment = labels[study][sample]['treatment'][index]
                if treatment is None:
                    if st:
                        break
                    continue
                if treatment == 'null':
                    labels[study][sample]['treatment'][index] = None
                    if st:
                        break
                    continue
                labels[study][sample]['treatment'][index] = labels[study][sample]['treatment'][index].replace(' ','_')
                labels[study][sample]['treatment'][index] = labels[study][sample]['treatment'][index].replace('_stress','')
                if 'nacl' in treatment:
                    labels[study][sample]['treatment'][index] = 'salt'
                if 'salt' in treatment:
                    labels[study][sample]['treatment'][index] = 'salt'
                if st:
                    break


            # if labels[study][sample]['tissue'] is None:
            #     labels[study][sample]['tissue'] = 'unknown'
            # for index in range(len(labels[study][sample]['tissue'])):
            #     labels[study][sample]['tissue'] = labels[study][sample]['tissue'].replace(' ','_')
            #     labels[study][sample]['tissue'] = labels[study][sample]['tissue'].replace('plants','plant')
            #     labels[study][sample]['tissue'] = labels[study][sample]['tissue'].replace('leaves','leaf')
            #     labels[study][sample]['tissue'] = labels[study][sample]['tissue'].replace('rossete','rosette')
            #     tissue = labels[study][sample]['tissue']
            #     if 'rosette' in tissue:
            #         labels[study][sample]['tissue'] = 'rosette'
                
            #     if 'seedling' in tissue:
            #         labels[study][sample]['tissue'] = 'seedling'
            # labels[study][sample]['tissue']
    return labels

def save_json(object,name, path):
    try:
        with open(f'{path}/{name}.json', 'w') as handle:
            json.dump(object, handle)
    except:
        pass

#cluster only the data of the most active genes

def make_df_from_labels(data_dict,col):
    records = []
# Iterate through the dictionary structure
    for study_id, samples in data_dict.items():
        for sample_id, sample_data in samples.items():
            # Create a new record with study_id included
            record = {
                'study_id': study_id,
                'sample_id': sample_id,
                **sample_data  # Unpack all sample data
            }
            records.append(record)

    # Create DataFrame
    df = pd.DataFrame(records)

    # Set sample_id as index
    df.set_index('sample_id', inplace=True)

    # Select and order the columns you want
    # df = df[['study_id', 'control', 'treatment', 'tissue', 'agar']]
    return df[col]

def keys_upper(test_dict):
    res = dict()
    for key in test_dict.keys():
        if isinstance(test_dict[key], dict):
            res[key.upper()] = keys_upper(test_dict[key])
        else:
            res[key.upper()] = test_dict[key]
    return res
def get_biggest_studies(labels,sample_index,k):
    studies_map,counts =np.unique(np.array(labels),return_counts = True)
    s = list(zip(studies_map,counts))
    s.sort(key=lambda x: x[1],reverse=True)
    s = s[:k]
    s= list(map(lambda x: x[0],s))
    mask =list(map(lambda x: x in s,labels))#TODO: make this a mask with true being the largest studies
    studies = np.unique(np.array(labels)[list(mask)])
    samples = []
    for s in studies:
        index = labels.index(s)
        samples.append((sample_index[index],s))
    x = 0
    return samples
def get_incomplete_studies(maps,sample_index,label):
    mask = list(map(lambda x: x is None,maps[label]))
    studies = np.unique(np.array(maps['study'])[list(mask)])
    samples = []
    for s in studies:
        index = maps['study'].index(s)
        samples.append(sample_index[index])
    x = 0
    return samples

def load_labels_study(path):
    labels = {}
    for file in os.listdir(path):
        if file.endswith('.json'):
            file_path = os.path.join(path, file)
            study = file_path.split('/')[-1].split('.')[0]
            if re.search('GSM',study):
                continue
            with open(file_path, 'r') as file:
                data = json.load(file)
            labels[study] = data
    return labels


def plot_pie_chart(data,path):
    sqrt = math.ceil(math.sqrt(len(data)))
    fig, axes = plt.subplots(sqrt, sqrt, figsize=(12, 10))
    # if sqrt*sqrt > len(data):
    #     to_remove = sqrt*sqrt- len(data)
    #     for i in range(to_remove):
    #         fig.delaxes(axes[sqrt-1,(sqrt-1)-i]) # The indexing is zero-based here
    axes = axes.ravel()
    for i, (key, values) in enumerate(data.items()):
        values = list(map(lambda x: str(x), values))
        unique, counts = np.unique(values, return_counts=True)
        axes[i].pie(counts, labels=unique, autopct='%1.1f%%', startangle=90)
        axes[i].set_title(f'Label Distribution - {key}')
        axes[i].axis('equal')

    plt.tight_layout()
    plt.savefig(f'{path}/pie.svg')
    plt.close()
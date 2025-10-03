import json
import os
from tqdm import tqdm
from condense_and_evaluate_labels import compare_labels
# from extraction_test_core import *
from meta_data_processing.utils.extractors_full import *
from llm_utils import get_condensed_labels
from classes import LabelMap
# TODO
#repeated function:

def load_labels_study(path):
    labels = {}
    for file in os.listdir(path):
        if file.endswith('.json'):
            file_path = os.path.join(path, file)
            study = file_path.split('/')[-1].split('.')[0]
            with open(file_path, 'r') as file:
                data = json.load(file)
            labels[study] = data
    return labels

def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object


def condense_labels(Studies,in_folder, saving_path,llm_grounding:bool = True):
    os.makedirs(saving_path,exist_ok=True)
    try:
        labels = load_labels_study(saving_path)
        raise ValueError()
    except:
        labels= {}
    

    # with open('constants/sample_list.json', 'r') as file:
    #     list_of_samples = json.load(file)
    seen = LabelMap('maps')
    count = 0
    for file in tqdm(os.listdir(in_folder)):
        if file.endswith('.json'):
            file_path = os.path.join(in_folder, file)
            
            
            study_id = file_path.split('/')[-1].split('.')[0]
            
            study_info = load_json(file_path)
            for sample_id in study_info:
                if sample_id == 'Study metadata':
                    continue
                if not (study_id in Studies):
                    continue
                # if  study_id in labels:
                #     if sample_id in labels[study_id]:
                #         continue
                # if not (sample_id in ['GSM2406811']):# ,'GSM2109731'
                    # continue
                
                sample_info = study_info[sample_id]
                # if study_id == 'GSE41935':
                #     temp.append(sample_id)
                fun = eval(f'{study_id}_extractor')
                python_object = fun(sample_info)
                # add the possibly missing control treatments:
                if python_object['treatment'] == []:
                    python_object['treatment'] = ['no treatment/control']
                x=0
                if llm_grounding and seen.check_past(python_object):
                    condensed = dict(get_condensed_labels(study_info=study_info['Study metadata'], sample_info=python_object))
                    seen.add_mapping(python_object,condensed)
                    python_object = condensed
                else:
                    condensed = seen.apply_mappings(python_object)
                if not (study_id in labels):
                    labels[study_id] = {}
                condensed['treatment'] = list(set(condensed['treatment']))
                labels[study_id][sample_id] = dict(condensed)
                seen.save_map()

    for study in labels:
        with open(f'{saving_path}/{study}.json', 'w') as handle:
            json.dump(labels[study], handle)


if __name__ == '__main__':
    Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']
    # Studies = ['GSE34188']
    in_folder = 'study_batch_metadata'
    experiment = 'venice'
    model = 'extractors_and_gemini'
    saving_path = f'labels/{model}/{experiment}'
    condense_labels(Studies,in_folder,saving_path)

    labels_1 = load_labels_study(saving_path)
    res = []
    for st in labels_1:
        for sam in labels_1[st]:
            save = labels_1[st][sam]
            save['id'] = sam
            res.append(labels_1[st][sam])
            
    with open(f'llm_condensed_labels.json', 'w') as handle:
        json.dump(res, handle)
    labels_2 = load_labels_study('/home/alex/Documents/GitHub/meta_data/in_use_labels/sample_labels_gemma_api_3_clean')
    compare_labels(labels_1,labels_2)

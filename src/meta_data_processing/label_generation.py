import json
import os
from tqdm import tqdm
# from extraction_test_core import *
import sys
module_dir = './'
sys.path.append(module_dir)
from src.meta_data_processing.utils.extractors_full import *
from src.meta_data_processing.utils.llm_utils import get_condensed_labels
from src.meta_data_processing.utils.classes import LabelMap
from src.constants import *
from src.meta_data_processing.utils.llm_utils import get_condensed_labels
from src.meta_data_processing.utils.classes import LabelMap


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


def condense_labels(Studies=Studies,in_folder=f'{METADATA_OUTPUT_DIR}/study_batch_metadata/', saving_path=LABELS_PATH,llm_grounding:bool = True):
    os.makedirs(saving_path,exist_ok=True)
    labels= {}
    seen = LabelMap('./data/maps')
    for file in tqdm(os.listdir(in_folder)):
        if file.endswith('.json'):
            file_path = os.path.join(in_folder, file)
            
            
            study_id = file_path.split('/')[-1].split('.')[0]
            
            study_info = load_json(file_path)
            for sample_id in study_info:
                if sample_id == 'study_metadata':
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
                if   "Light 24h" in python_object['treatment']:
                    print(f'study: {study_id} sample: {sample_id}')
                if llm_grounding and seen.check_past(python_object):
                    condensed = dict(get_condensed_labels(study_info=study_info['study_metadata'], sample_info=python_object))
                    seen.add_mapping(python_object,condensed)
                    python_object = condensed
                else:
                    condensed = seen.apply_mappings(python_object)
                if not (study_id in labels):
                    labels[study_id] = {}
                def flatten(lst,ret = []):
                    for el in lst:
                        if type(el) == list:
                            flatten(el,ret)
                        else:
                            ret.append(el)
                    return ret

                condensed['treatment'] = list(set(flatten(condensed['treatment'])))
                labels[study_id][sample_id] = dict(condensed)
                seen.save_map()

    for study in labels:
        with open(f'{saving_path}/{study}.json', 'w') as handle:
            json.dump(labels[study], handle)


if __name__ == '__main__':
    condense_labels()

    labels_1 = load_labels_study(LABELS_PATH)
    res = []
    for st in labels_1:
        for sam in labels_1[st]:
            save = labels_1[st][sam]
            save['id'] = sam
            res.append(labels_1[st][sam])
            
    with open(f'llm_condensed_labels.json', 'w') as handle:
        json.dump(res, handle)

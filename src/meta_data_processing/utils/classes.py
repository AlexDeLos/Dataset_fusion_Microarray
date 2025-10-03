from typing import List, Optional
import json
import copy

def load_json(path:str):
    with open(path, 'r') as file:
        object = json.load(file)
    return object


class LabelMap:
    def __init__(self, path:Optional[str]=None):
        self.path = path
        if path is None:
            self.map_treatment = {}
            self.map_tissue = {}
            self.map = {}
        else:
            try:
                self.map_treatment = load_json(path+'/map_treatment.json')
                self.map_tissue = load_json(path+'/map_tissue.json')
                self.map = load_json(path+'/map.json')
            except:
                Warning('Path not found')
                self.map_treatment = {}
                self.map_tissue = {}
                self.map = {}

    def add(self,label:str,id)->None:
        self.map[label] = id

    def add_treatment(self,label:str,id)->None:
        self.map_treatment[label] = id

    def add_tissue(self,label:str,id)->None:
        self.map_tissue[label] = id
    
    def add_mapping(self,og,grounded)->None:
        for el in og:
            if  isinstance(og[el],List):
                for i,term in enumerate(og[el]):
                    if len(og[el]) != len(grounded[el]):
                        self.add_treatment(str(og[el]),grounded[el])
                        break
                    else:
                        self.add_treatment(term,grounded[el][i])
                        
            elif (el =='tissue'):
                self.add_tissue(og[el],grounded[el])
            else:
                self.add(og[el],grounded[el])

    def in_maps_evaluated(self,el)->bool:
        return el in self.map_treatment or el in self.map_tissue

    def in_maps(self,el)->bool:
        return self.in_maps_evaluated(el) or el in self.map
    

    def check_past(self,sample:dict)->bool:
        run = False
        for el in sample:
            if  isinstance(sample[el],List):
                for term in sample[el]:
                    if not self.in_maps(term):
                        if self.in_maps(str(sample[el])):
                            pass
                        else:
                            run = True
            elif el !='id':
                if not self.in_maps(sample[el]):
                    run = True
        return run
    
    def save_map(self) -> None:
        with open(f'{self.path}/map_treatment.json', 'w') as handle:
            json.dump(self.map_treatment, handle)
        with open(f'{self.path}/map_tissue.json', 'w') as handle:
            json.dump(self.map_tissue, handle)
        with open(f'{self.path}/map.json', 'w') as handle:
            json.dump(self.map, handle)
    
    
    def apply_mappings(self,og:dict)-> dict:
        ret = copy.deepcopy(og)
        for el in og:
            if  isinstance(og[el],List):
                for i,_ in enumerate(og[el]):
                    try:
                        ret[el][i] = self.map_treatment[og[el][i]]
                    except:
                        ret[el] =self.map_treatment[str(og[el])]
                        break
   
            elif el =='tissue':
                    ret[el] = self.map_tissue[og[el]]
            else:
                    ret[el] = self.map[og[el]]

        return ret
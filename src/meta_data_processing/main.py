# Gather the metadata
from Generate_sample_metadata import generate_metadata

# Studies = ['GSE44053', 'GSE77815', 'GSE16474', 'GSE4062', 'GSE9415', 'GSE18624', 'GSE40061', 'GSE112161', 'GSE20494', 'GSE79997', 'GSE27552', 'GSE16222', 'GSE110857', 'GSE4760', 'GSE46205', 'GSE71001', 'GSE58616', 'GSE162310', 'GSE22107', 'GSE62163', 'GSE51897', 'GSE72949', 'GSE90562', 'GSE5628', 'GSE26266', 'GSE34188', 'GSE34595', 'GSE76827', 'GSE119383', 'GSE65046', 'GSE11758', 'GSE65414', 'GSE37408', 'GSE5624', 'GSE10643', 'GSE15577', 'GSE11538', 'GSE70861', 'GSE6583', 'GSE27551', 'GSE12619', 'GSE121225', 'GSE108070', 'GSE78713', 'GSE110079', 'GSE63128', 'GSE60960', 'GSE37118', 'GSE79681', 'GSE63372', 'GSE5622', 'GSE26983', 'GSE27550', 'GSE19603', 'GSE95202', 'GSE53308', 'GSE16765', 'GSE71855', 'GSE58620', 'GSE24177', 'GSE35258', 'GSE10670', 'GSE49418', 'GSE18666', 'GSE83136', 'GSE44655', 'GSE27549', 'GSE19700', 'GSE103398', 'GSE63522', 'GSE201609', 'GSE5620', 'GSE66369', 'GSE2268', 'GSE71237', 'GSE48474', 'GSE41935', 'GSE27548', 'GSE5623', 'GSE72050', 'GSE126373']
Studies = ['GSE201609','GSE46205','GSE27548','GSE41935','GSE34188','GSE27550','GSE76827','GSE5620','GSE110079','GSE60960','GSE5624','GSE5622','GSE126373','GSE5628','GSE63128']
# generate_metadata()


# Condese the metadata into the give sample labels
from label_generation import condense_labels, load_labels_study
model = 'gemini-2.5-flash'

# Studies = ['GSE34188']
in_folder = 'sample_metadata'
experiment = 'total_Sumer_school'

saving_path = f'labels/{model}/{experiment}'

condense_labels(Studies=Studies, in_folder=in_folder,model=model,saving_path=saving_path)
# Check the label diference with the other sets of labels

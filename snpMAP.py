# This will be the startup script
import os

import numpy as np

from utils import parser
from utils.dir_utils import prepare_directiories, clear_cache
from features.feature_extactor import FeatureExtractor
from preprocessing import preprocess
from classification.model import SnpMAPModel, model_trained
from io_utils.writing import ResultsWriter
from constants import default_mv_fields


config = parser.load_config()
cache_dir = config['output']['cache_dir']
results_dir = config['output']['results_dir']

fpaths = {
    'train': config['input']['clinvar_path'],  # does not need to be shrunken or annotated
    'work': config['input']['case_path']
}


def run(mode, SNP_type, classification):

    shrunken_fpath = preprocess.shrink_vcf(fpaths[mode], mode)

    annotated_fpath = preprocess.annotate_vcf(shrunken_fpath, mode)

    SNPs_path = preprocess.extract_SNPs(
        annotated_fpath, SNP_type=SNP_type, mode=mode)

    # hgvs_fpath = preprocess.extract_hgvs_ids(SNPs_path)

    labels_fpath = preprocess.extract_labels(SNPs_path, mode, classification)

    fe = FeatureExtractor(
        SNPs_path, SNP_type=SNP_type, mode=mode,
        mv_fields=default_mv_fields[SNP_type])

    fe.construct_features()

    features_fpath, features_labels_fpath = fe.make_dataset(
        labels_fpath, classification)

    # print(features_fpath, features_labels_fpath)
    model = SnpMAPModel(features_fpath, features_labels_fpath, mode=mode, SNP_type=SNP_type)
        
    if mode == 'work':

        results = model.predict()

        return results['results']

    elif mode == 'train':

        model.train()

        return None


def main():
    
    prepare_directiories(cache_dir, results_dir)
    results = []
    for SNP_type in ['nsSNP', 'sSNP']: 
        classification = config['classification']
        if config['retrain'] == 'True' or not model_trained(SNP_type):
            print('Preprating to train the model.')
            # clear_cache('train', cache_dir)
            run('train', SNP_type, classification)
            print('Model training DONE.')

        print(f'Preparing to classify the {SNP_type} variants the case file.')
        r = run('work', SNP_type, classification)
        print(f'Prediction for the type {SNP_type} finished.')
        results.append(r)

    print(np.unique(list(results[0].values()), return_counts=True))
    print(np.unique(list(results[1].values()), return_counts=True))

    # print(len(results[0]))
    # print(len(results[1]))
    R = results[1]
    R.update(results[0])
    # print(len(R))
    writer = ResultsWriter(fpaths['work'])
    writer.write_results(R)

    print('DONE')


if __name__ == "__main__":
    main()

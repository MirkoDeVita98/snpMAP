from time import time
import os.path
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import stats
from Bio.Data.IUPACData import protein_letters_3to1
import myvariant
import rhapsody as rd

from utils import parser, file_info
from constants import default_mv_fields, NS_SNP, S_SNP, I_SNP
from preprocessing.preprocess import get_hgvs_ids

from dependencies.seattle138.query_seattle import download_seattle


import pickle


config = parser.load_config()
cache_dir = config['output']['cache_dir']
use_polyphen = config['use_polyphen']


class FeatureExtractor:
    def __init__(self, vcf_fpath, mode, SNP_type, mv_fields):
        self.fields = mv_fields
        self.vcf_fpath = vcf_fpath
        self.hgvs_ids = get_hgvs_ids(self.vcf_fpath)
        self.mv = myvariant.MyVariantInfo()
        self.mv_results = None
        self.mode = mode
        self.SNP_type = SNP_type
        self.uniprot_ids_fpath = os.path.join(
            cache_dir, f'uniprot_{self.SNP_type}_{self.mode}.csv')
        self.uniprot_ids_mapper_fpath = os.path.join(
            cache_dir, f'uniprot_mapper_{self.SNP_type}_{self.mode}.csv')
        if SNP_type == NS_SNP:
            self.polyphen2_fpath = os.path.join(cache_dir,
                f'polyphen2_{self.SNP_type}_{self.mode}.csv')
                # f'polyphen2_{self.mode}.csv')
            self.BLOSUM_fpath = os.path.join(cache_dir,
                f'BLOSUM_{self.SNP_type}_{self.mode}.csv')
            self.conservation_fpath = os.path.join(cache_dir,
                f'conservation_{self.SNP_type}_{self.mode}.csv')
            self.gerp_fpath = os.path.join(cache_dir,
                f'gerp_{self.SNP_type}_{self.mode}.csv')
            self.sig_fpath = os.path.join(cache_dir,
                f'sig_{self.SNP_type}_{self.mode}.csv')
            self.gc_fpath = os.path.join(cache_dir,
                f'gc_{self.SNP_type}_{self.mode}.csv')
            self.encode_fpath = os.path.join(cache_dir,
                f'encode_{self.SNP_type}_{self.mode}.csv')
            self.dist_fpath = os.path.join(cache_dir,
                f'dist_{self.SNP_type}_{self.mode}.csv')
        elif SNP_type == S_SNP:
            self.seattle_download_fpath = os.path.join(cache_dir,
                f'seattle_annotation_{self.SNP_type}_{self.mode}.vcf')
            self.conservation_fpath = os.path.join(cache_dir,
                f'conservation_{self.SNP_type}_{self.mode}.csv')

            self.gerp_fpath = os.path.join(cache_dir,
                f'gerp_{self.SNP_type}_{self.mode}.csv')

            self.encode_fpath = os.path.join(cache_dir,
                f'encode_{self.SNP_type}_{self.mode}.csv')

            self.gc_fpath = os.path.join(cache_dir,
                f'gc_{self.SNP_type}_{self.mode}.csv')

            self.cpg_fpath = os.path.join(cache_dir,
                f'cpg_{self.SNP_type}_{self.mode}.csv')

            self.dist_fpath = os.path.join(cache_dir,
                f'dist_{self.SNP_type}_{self.mode}.csv')

            self.seattle_fpath = os.path.join(cache_dir,
                f'seattle_{self.SNP_type}_{self.mode}.csv')
        self.uniprot2hgvs = None

    
    def construct_features(self):
        if self.SNP_type == NS_SNP:
            self.download_features()

            self.extract_uniprot_ids()

            self.extract_conservation_features()
            if config['use_polyphen'] == 'True':
                self.extract_polyphen2_features()
            self.extract_BLOSUM_scores()
            self.extract_gc_features()
            self.extract_encode_features()
            self.extract_dist_features()
            self.extract_gerp_features()
            self.extract_sig_features()
        elif self.SNP_type == S_SNP:
            self.download_features()
            self.download_seattle_features()
            self.extract_conservation_features()
            self.extract_gerp_features()
            self.extract_dist_features()
            self.extract_gc_features()
            self.extract_cpg_features()
            self.extract_encode_features()
            self.extract_seattle_features()

    def download_features(self):

        results_fpath = os.path.join(
            cache_dir, f'mv_results_{self.SNP_type}_{self.mode}.pickle')

        if self.mv_results is not None:
            print('MyVariant features already downloaded.')
            return

        if os.path.exists(results_fpath):
            print('MyVariant features already downloaded.')
            self.mv_results = pickle.load(open(results_fpath, 'rb'))
            return

        print(f'Downloading the following features from MyVariant: {self.fields}')
        print(f'Number of SNPs to query for: {len(self.hgvs_ids)}.')
        start_time = time()
        self.mv_results = self.mv.getvariants(
            self.hgvs_ids, fields=self.fields)
        end_time = time()
        elapsed_time = end_time - start_time
        print(f'\nMyVariant feature download finished in {int(elapsed_time // 60)} minutes and {int(elapsed_time % 60)} seconds.')

        pickle.dump(self.mv_results, open(results_fpath, 'wb'))


    def extract_uniprot_ids(self):
        # if os.path.exists(self.uniprot_ids_fpath):
        #     return

        print('Extracting UniProt ids from MyVariant results.')
        uniprot_lines = []
        uniprot_dict = {}
        for ii, result in tqdm(
            enumerate(self.mv_results), total=len(self.mv_results)):
            lines = self._get_uniprot_lines(result)
            if len(lines) > 0:
                for line in lines:
                    uniprot_dict[line] = self.hgvs_ids[ii]
                uniprot_lines.extend(lines)
        
        with open(self.uniprot_ids_fpath, 'w') as fw:
            fw.write('\n'.join(uniprot_lines) + '\n')

        self.uniprot2hgvs = uniprot_dict
        print('Extracting UniProt ids from MyVariant results finished.')
        
    
    def extract_gc_features(self):
        if os.path.exists(self.gc_fpath):
            return
        F = []
        for result in self.mv_results:
            f = self._get_gc_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(F, columns=['ID', 'gc'])
        df = df.set_index(['ID'])
        df.to_csv(self.gc_fpath)
    
    def extract_cpg_features(self):
        if os.path.exists(self.cpg_fpath):
            return
        F = []
        for result in self.mv_results:
            f = self._get_cpg_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(F, columns=['ID', 'cpg'])
        df = df.set_index(['ID'])
        df.to_csv(self.cpg_fpath)

    def extract_dist_features(self):
        if os.path.exists(self.dist_fpath):
            return
        print('Extracting distance to closest transcribed sequence.')
        F = []
        for result in tqdm(self.mv_results):
            f = self._get_dist_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(F, columns=['ID', 'min_dist_tss', 'min_dist_tse'])
        df = df.set_index(['ID'])
        df.to_csv(self.dist_fpath)
        print('Extracting distance to closest transcribed sequence features finished.')


    def extract_encode_features(self):
        if os.path.exists(self.encode_fpath):
            return
        print('Extracting encoding features.')
        F = []
        for result in tqdm(self.mv_results):
            f = self._get_encode_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(F, columns=[
            'ID', 'exp', 'h3k27ac', 'h3k4me1', 'h3k4me3', 'nucleo'])
        df = df.set_index(['ID'])
        df.to_csv(self.encode_fpath)
        print('Extracting encoding features finished.')


    def extract_gerp_features(self):
        if os.path.exists(self.gerp_fpath):
            return
        print('Extracting gerp features.')
        F = []
        for result in tqdm(self.mv_results):
            f = self._get_gerp_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(F, columns=['ID', 'n', 'rs', 's'])
        df = df.set_index(['ID'])
        df.to_csv(self.gerp_fpath)
        print('Extracting gerp features finished.')
    

    def extract_sig_features(self):
        if os.path.exists(self.sig_fpath):
            return
        print('Extracting sig features.')
        F = []
        for result in tqdm(self.mv_results):
            f = self._get_sig_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(
            F, columns=['ID', 'ctcf', 'dnase', 'faire', 'myc', 'polii'])
        df = df.set_index(['ID'])
        df.to_csv(self.sig_fpath)
        print('Extracting sig features finished.')


    def extract_conservation_features(self):

    
        if os.path.exists(self.conservation_fpath):
            return
        print('Extracting conservation features.')
        F = []
        for result in tqdm(self.mv_results):
            f = self._get_conservation_features(result)
            if len(f) > 0:
                F.append(f)

        df = pd.DataFrame(
            F, columns=[
                'ID', 'phast_cons_mammalian', 'phast_cons_primate',
                'phast_cons_vertebrate', 'phylop_mammalian', 'phylop_primate',
                'phylop_vertebrate'])
        df = df.set_index(['ID'])
        df.to_csv(self.conservation_fpath)
        print('Extracting conservation features finished.')

    def download_seattle_features(self):
        if os.path.exists(self.seattle_download_fpath):
            print('SeattleSeq features already downloaded.')
            return
        
        print("Downloading SeattleSeq features...")
        download_seattle(config["SeattleSeq"]["java_classpath"],self.vcf_fpath,self.seattle_download_fpath)


    def extract_seattle_features(self):
        num_lines = sum(1 for line in open(self.seattle_download_fpath,'r'))
        with open(self.seattle_download_fpath,'r') as f, open(self.seattle_fpath,'w') as fw:
            fw.write("ID,d2splice,tfbs\n")
            last_processed = ""
            for line in tqdm(f, total=num_lines):
                if(line[0] != '#'):
                    x= line.split('\t')
                    id = 'chr'+ x[1] +':g.' + x[2] + x[3] +'>'  
                    if x[5][0] == x[3]:
                        id += x[5][2]
                    else:
                        id += x[5][0] 
                    d2splice = x[12]
                    if  x[13] == "none" or  x[13] == "none\n":
                        tfbs = 0
                    else:
                        tfbs = 1
                    if id != last_processed:
                        last_processed = id
                        fw.write(id + ","+ str(d2splice) + "," + str(tfbs) +"\n")


    def download_polyphen_features(self):

        if not use_polyphen == 'True':
            print('Will not download Polyphen features.')
            return

        if os.path.exists(os.path.join(
            cache_dir, f'pph2-res-{self.SNP_type}-{self.mode}-full.txt')):
            print('Polyphen results have alrady been downloaded.')
            return

        print('Downloading PolyPhen2 features.')
        print('Depending on the server load, this might take some time.')
        print('If this process hangs, it may work later, or you can set the ' +
              '"use_polyphen2" flag to False in the config file to skip this.')

        try:
            rd.features.PolyPhen2.queryPolyPhen2(
                self.uniprot_ids_fpath,
                prefix=os.path.join(cache_dir,
                    f'pph2-res-{self.SNP_type}-{self.mode}'),
                ignore_errors=True)
        except:
            pass


    def extract_polyphen2_features(self):

        if os.path.exists(self.polyphen2_fpath):
            print('Polyphen features havel already been extracted.')
            return

        self.download_polyphen_features()

        pp2_results = pd.read_csv(
            os.path.join(cache_dir,
            f'pph2-res-{self.SNP_type}-{self.mode}-full.txt'),
            sep='\t').dropna()
        pp2_results['   pos'] = pp2_results['   pos'].astype(np.int)
        pp2_results = pp2_results.set_index(
            ['acc       ', '   pos', 'aa1', 'aa2'])

        dataset = {
            'score1': {},
            'dscore': {},
            'nobs': {},
            # 'dvol': {},
            'transv': {},
            'cpg': {},
            'idpmax': {},
            # 'idpsnp': {},
            'idqmin': {},
            'mindjxn': {},
        }

        key_occurences = {
            'score1': defaultdict(int),
            'dscore': defaultdict(int),
            'nobs': defaultdict(int),
            # 'dvol': defaultdict(int),
            'transv': defaultdict(int),
            'cpg': defaultdict(int),
            'idpmax': defaultdict(int),
            # 'idpsnp': defaultdict(int),
            'idqmin': defaultdict(int),
            'mindjxn': defaultdict(int),
        }

        keys = set()
        for row in tqdm(pp2_results.iterrows(), total=pp2_results.shape[0]):
            key = self._get_key(*row[0])
            keys.add(key)

            self._parse_row_for_info(
                row[1], key_occurences['score1'], dataset['score1'], key,
                'Score1')
            self._parse_row_for_info(
                row[1], key_occurences['dscore'], dataset['dscore'], key,
                'dScore')
            self._parse_row_for_info(
                row[1], key_occurences['nobs'], dataset['nobs'], key,
                '  Nobs')
            self._parse_row_for_info(
                row[1], key_occurences['transv'], dataset['transv'], key,
                'Transv', dtype=np.int)
            self._parse_row_for_info(
                row[1], key_occurences['cpg'], dataset['cpg'], key,
                'CpG', dtype=np.int)
            self._parse_row_for_info(
                row[1], key_occurences['idpmax'], dataset['idpmax'], key,
                '  IdPmax')
            self._parse_row_for_info(
                row[1], key_occurences['idqmin'], dataset['idqmin'], key,
                '  IdQmin')
            self._parse_row_for_info(
                row[1], key_occurences['mindjxn'], dataset['mindjxn'], key,
                ' MinDJxn')

        common_keys = set()

        interesting_keys = [
            'score1',
            'dscore',
            'nobs',
            'transv',
            'cpg',
            'idpmax',
            'idqmin',
            'mindjxn'
        ]

        # for key in key_occurences['score1']:
        #     if all([key in key_occurences[d] for d in interesting_keys]):
        #         common_keys.add(key)

        # lP = []
        # for key in tqdm(common_keys, total=len(common_keys)):
        #     if key in self.uniprot2hgvs.keys():
        #         lP.append(
        #             np.asarray([
        #                 self.uniprot2hgvs[key],
        #                 dataset['score1'][key],
        #                 dataset['dscore'][key],
        #                 dataset['nobs'][key],
        #                 dataset['transv'][key],
        #                 dataset['cpg'][key],
        #                 dataset['idpmax'][key],
        #                 dataset['idqmin'][key],
        #                 dataset['mindjxn'][key]]))

        lP = []
        for key in keys:
            if key not in self.uniprot2hgvs:
                continue
            f = [self.uniprot2hgvs[key]]
            for feature in interesting_keys:
                f.append(dataset[feature][key]
                    if key in dataset[feature] else np.nan)
            lP.append(np.asarray(f))

        df = pd.DataFrame(
            np.asarray(lP),
            columns=[
                'ID', 'Score1', 'dScore', 'Nobs', 'Transv',
                'CpG', 'IdPmax', 'IdQmin', 'MinDJxn'])
        df = df.set_index(['ID'])

        df.to_csv(self.polyphen2_fpath)


    def extract_BLOSUM_scores(self):
        if os.path.exists(self.BLOSUM_fpath):
            return
        SAVs = np.unique(pd.read_csv(self.uniprot_ids_fpath, header=None)[0].to_list())
        bs = rd.features.BLOSUM.calcBLOSUMfeatures(SAV_coords=SAVs)
        hvgss = [self.uniprot2hgvs[uniprot] for uniprot in SAVs]
        df = pd.DataFrame().from_dict(
            {'ID': hvgss, 'blosum': [np.int(b[0]) for b in bs]})
        df = df.set_index('ID')
        df.to_csv(self.BLOSUM_fpath)


    def make_dataset(self, labels_fpath, classification):
        
        print('Making the final dataset.')

        if self.SNP_type == NS_SNP:
            res = self.make_dataset_nsSNP(labels_fpath, classification)
        if self.SNP_type == S_SNP:
            res = self.make_dataset_sSNP(labels_fpath, classification)

        print('Making the final dataset finished.')

        return res

    
    def make_dataset_sSNP(self, labels_fpath, classification):
        
        features_result_fpath = os.path.join(
                cache_dir, f'features_{self.SNP_type}_{self.mode}.csv')

        labels_result_fpath = None
        labels_result_fpath = os.path.join(
            cache_dir, f'features_labels_{self.SNP_type}_{self.mode}_{classification}.csv')

        if os.path.exists(features_result_fpath):
            print('The dataset has alrady been constructed.')
            return features_result_fpath, labels_result_fpath

        C = pd.read_csv(self.conservation_fpath).set_index(['ID'])
        G = pd.read_csv(self.gc_fpath).set_index(['ID'])
        D = pd.read_csv(self.dist_fpath).set_index(['ID'])
        E = pd.read_csv(self.encode_fpath).set_index(['ID'])
        Gr = pd.read_csv(self.gerp_fpath).set_index(['ID'])
        CPG = pd.read_csv(self.cpg_fpath).set_index(['ID'])
        S = pd.read_csv(self.seattle_fpath).set_index(['ID'])

        Ds = [ C, G, D, E, Gr,CPG,S]

        M = Ds[np.argmin(np.asarray([El.shape[0] for El in Ds]))]

        print([El.shape for El in Ds])

        if self.mode == 'train':
            labels = pd.read_csv(labels_fpath).set_index(['ID'])

        X = []
        y = []
        for ii, ix in tqdm(enumerate(M.index), total=M.shape[0]):
            if all([ix in El.index for El in Ds]):
                f = [ix]
                if len(C.loc[ix].shape) > 1:
                    f.extend(C.loc[ix].mean().to_numpy())
                else:
                    f.extend(C.loc[ix].to_numpy())

                if len(G.loc[ix].shape) > 1:
                    f.extend(G.loc[ix].mean().to_numpy())
                else:
                    f.extend(G.loc[ix].to_numpy())

                if len(D.loc[ix].shape) > 1:
                    f.extend(D.loc[ix].mean().to_numpy())
                else:
                    f.extend(D.loc[ix].to_numpy())

                if len(E.loc[ix].shape) > 1:
                    f.extend(E.loc[ix].mean().to_numpy())
                else:
                    f.extend(E.loc[ix].to_numpy())

                if len(Gr.loc[ix].shape) > 1:
                    f.extend(Gr.loc[ix].mean().to_numpy())
                else:
                    f.extend(Gr.loc[ix].to_numpy())
                
                if len(CPG.loc[ix].shape) > 1:
                    f.extend(CPG.loc[ix].mean().to_numpy())
                else:
                    f.extend(CPG.loc[ix].to_numpy())

                if len(S.loc[ix].shape) > 1:
                    f.extend(S.loc[ix].mean().to_numpy())
                else:
                    f.extend(S.loc[ix].to_numpy())

                if self.mode == 'train':
                    y.append((ix, labels.loc[ix][
                        0 if classification == 'binary' else 1])) # TODO

                X.append(f)
        
        cols = ['ID']
        for E in Ds:
            cols.extend(E.columns)
        df = pd.DataFrame(np.asarray(X), columns=cols)
        df = df.set_index('ID')
        df.to_csv(features_result_fpath)
        
        if self.mode == 'train':
            df = pd.DataFrame.from_dict(
                {'ID': [j[0] for j in y], 'label': [j[1] for j in y]})
            df = df.set_index('ID')
            df.to_csv(labels_result_fpath)

        return features_result_fpath, labels_result_fpath


    def make_dataset_nsSNP(self, labels_fpath, classification):

        if self.SNP_type == NS_SNP:
            features_result_fpath = os.path.join(
                cache_dir, f'features_{self.SNP_type}_{"polyphen" if use_polyphen == "True" else "no_polyphen"}_{self.mode}.csv')
        else:
            features_result_fpath = os.path.join(
                cache_dir, f'features_{self.SNP_type}_{self.mode}.csv')

        binary_labels_result_fpath = None
        multiclass_labels_result_fpath = None
        if self.mode == 'train':
            if self.SNP_type == NS_SNP:
                binary_labels_result_fpath = os.path.join(
                    cache_dir, f'features_labels_{self.SNP_type}_{"polyphen" if use_polyphen == "True" else "no_polyphen"}_{self.mode}_binary.csv')
                multiclass_labels_result_fpath = os.path.join(
                    cache_dir, f'features_labels_{self.SNP_type}_{"polyphen" if use_polyphen == "True" else "no_polyphen"}_{self.mode}_multiclass.csv')
                
            else:
                binary_labels_result_fpath = os.path.join(
                    cache_dir, f'features_labels_{self.SNP_type}_{self.mode}_binary.csv')
                multiclass_labels_result_fpath = os.path.join(
                    cache_dir, f'features_labels_{self.SNP_type}_{self.mode}_multiclass.csv')
            
        labels_result_fpath = binary_labels_result_fpath if classification == 'binary' else multiclass_labels_result_fpath
        
        if os.path.exists(features_result_fpath):
            print('The dataset has alrady been constructed.')
            return features_result_fpath, labels_result_fpath

        B = pd.read_csv(self.BLOSUM_fpath).set_index(['ID'])
        C = pd.read_csv(self.conservation_fpath).set_index(['ID'])
        if config['use_polyphen'] == 'True':
            P = pd.read_csv(self.polyphen2_fpath).set_index(['ID'])
        G = pd.read_csv(self.gc_fpath).set_index(['ID'])
        D = pd.read_csv(self.dist_fpath).set_index(['ID'])
        E = pd.read_csv(self.encode_fpath).set_index(['ID'])
        Gr = pd.read_csv(self.gerp_fpath).set_index(['ID'])
        # S = pd.read_csv(self.sig_fpath).set_index(['ID'])
        # Ds = [B, C, P, G, D, E, Gr, S]
        if config['use_polyphen'] == 'True':
            Ds = [B, C, P, G, D, E, Gr]
        else:
            Ds = [B, C, G, D, E, Gr]

        M = Ds[np.argmin(np.asarray([El.shape[0] for El in Ds]))]

        print([El.shape for El in Ds])

        if self.mode == 'train':
            labels = pd.read_csv(labels_fpath).set_index(['ID'])

        X = []
        yb, ym = [], []
        for ii, ix in tqdm(enumerate(M.index), total=M.shape[0]):
            if all([ix in El.index for El in Ds]):
                f = [ix]

                if config['use_polyphen'] == 'True':
                    if len(P.loc[ix].shape) > 1:
                        f.extend(P.loc[ix].mean().to_numpy())
                    else:
                        f.extend(P.loc[ix].to_numpy())

                if len(B.loc[ix].shape) > 1:
                    # if len(np.unique(B.loc[ix]['blosum'])) != 1:
                    #     print('WARNING')
                    # f.append(B.loc[ix]['blosum'][0])  # TODO
                    f.append(stats.mode(B.loc[ix]['blosum'])[0][0])  # TODO
                else:
                    f.append(B.loc[ix]['blosum'])

                if len(C.loc[ix].shape) > 1:
                    f.extend(C.loc[ix].mean().to_numpy())
                else:
                    f.extend(C.loc[ix].to_numpy())

                if len(G.loc[ix].shape) > 1:
                    f.extend(G.loc[ix].mean().to_numpy())
                else:
                    f.extend(G.loc[ix].to_numpy())

                if len(D.loc[ix].shape) > 1:
                    f.extend(D.loc[ix].mean().to_numpy())
                else:
                    f.extend(D.loc[ix].to_numpy())

                if len(E.loc[ix].shape) > 1:
                    f.extend(E.loc[ix].mean().to_numpy())
                else:
                    f.extend(E.loc[ix].to_numpy())

                if len(Gr.loc[ix].shape) > 1:
                    f.extend(Gr.loc[ix].mean().to_numpy())
                else:
                    f.extend(Gr.loc[ix].to_numpy())

                # if len(S.loc[ix].shape) > 1:
                #     f.extend(S.loc[ix].mean().to_numpy())
                # else:
                #     f.extend(S.loc[ix].to_numpy())

                if self.mode == 'train':
                    yb.append((ix, labels.loc[ix][0])) # TODO
                    ym.append((ix, labels.loc[ix][1])) # TODO

                X.append(f)
        
        cols = ['ID']
        for E in Ds:
            cols.extend(E.columns)
        df = pd.DataFrame(np.asarray(X), columns=cols)
        df = df.set_index('ID')
        df.to_csv(features_result_fpath)
        
        if self.mode == 'train':
            df = pd.DataFrame.from_dict(
                {'ID': [j[0] for j in yb], 'label': [j[1] for j in yb]})
            df = df.set_index('ID')
            df.to_csv(binary_labels_result_fpath)

            df = pd.DataFrame.from_dict(
                {'ID': [j[0] for j in ym], 'label': [j[1] for j in ym]})
            df = df.set_index('ID')
            df.to_csv(multiclass_labels_result_fpath)
            
        return features_result_fpath, labels_result_fpath


    def make_labels_sSNP(self, labels_fpath):
        G = pd.read_csv(self.cadd_fpath).set_index(['ID'])
        C = pd.read_csv(self.conservation_fpath).set_index(['ID'])

        labels = pd.read_csv(labels_fpath).set_index(['ID'])

        y = []
        for ii, ix in tqdm(enumerate(G.index), total=G.shape[0]):
            if ix in G.index  and ix in C.index and ix in labels.label:
                y.append(labels.loc[ix].label)
        
        df = pd.DataFrame(np.asarray(y))
        df.index.name = 'ix'
        result_fpath = os.path.join(cache_dir, f'features_labels_{self.mode}.csv')
        df.to_csv(result_fpath)
        return result_fpath


    def _parse_row_for_info(self, row, keys_dict, val_dict, row_key,
                            column_key, dtype=np.float):

        if '?' not in str(row[column_key]):
            if keys_dict[row_key] > 0:
                o = val_dict[row_key] * keys_dict[row_key] \
                    if val_dict[row_key] != np.nan else 0
                n = (o + dtype(row[column_key])) / (keys_dict[row_key] + 1)

                val_dict[row_key] = n
            else:
                val_dict[row_key] = dtype(row[column_key])
            keys_dict[row_key] += 1


    def _get_key(self, pid, pos, a1, a2):
        pid = ''.join(pid.split())
        a1 = ''.join(a1.split())
        a2 = ''.join(a2.split())
        return f'{pid} {pos} {a1} {a2}'


    def _get_gerp_features(self, result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']

        if 'gerp' not in cadd:
            return []
        
        gerp = cadd['gerp']
        # print(gerp)
        # if any([key not in gerp for key in ['n', 'rs', 's']]):
        #     return []

        features = [result['_id']]
        keys = ['n', 'rs', 's']
        for key in keys:
            features.append(gerp[key] if key in gerp else np.nan)
        # features = [result['_id'], gerp['n'], gerp['rs'], gerp['s']]

        return features


    def _get_sig_features(self, result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']

        if 'encode' not in cadd or 'sig' not in cadd['encode']:
            return []
        
        sig = cadd['encode']['sig']

        keys = ['ctcf', 'dnase', 'faire', 'myc', 'polii']
        # if any([key not in sig for key in keys]):
        #     return []

        features = [result['_id']]
        for key in keys:
            features.append(sig[key] if key in sig else np.nan)
        # features.extend([sig[key] for key in sig])

        return features


    def _get_gc_features(self, result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']

        if 'gc' not in cadd:
            return []
        
        features = [result['_id'], cadd['gc']]

        return features

    def _get_cpg_features(self,result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']

        if 'cpg' not in cadd:
            return []
        
        features = [result['_id'], cadd['cpg']]

        return features

    def _get_dist_features(self, result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']

        # if 'min_dist_tss' not in cadd or 'min_dist_tse' not in cadd:
        #     return []
        
        features = [
            result['_id'], 
            cadd['min_dist_tss'] if 'min_dist_tss' in cadd else np.nan,
            cadd['min_dist_tse'] if 'min_dist_tse' in cadd else np.nan
        ]

        return features


    def _get_encode_features(self, result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']

        if 'encode' not in cadd:
            return []
        
        encode = cadd['encode']

        # if any([key not in encode for key in
        #     ['exp', 'h3k27ac', 'h3k4me1', 'h3k4me3', 'nucleo']]):
        #     return []
        keys = ['exp', 'h3k27ac', 'h3k4me1', 'h3k4me3', 'nucleo']
        features = [result['_id']]
        for key in keys:
            features.append(encode[key] if key in encode else np.nan)

        # features = [
        #     result['_id'], encode['exp'], encode['h3k27ac'], 
        #     encode['h3k4me1'], encode['h3k4me3'], encode['nucleo']]

        return features


    def _get_conservation_features(self, result):
        if 'cadd' not in result:
            return []

        cadd = result['cadd']
        # print(cadd)

        if 'phylop' not in cadd or 'phast_cons' not in cadd:
            return []
        
        phast_cons = cadd['phast_cons']
        phylop = cadd['phylop']

        keys = ['mammalian', 'primate', 'vertebrate']
        if any([key not in phast_cons for key in keys]):
            return []
        if any([key not in phylop for key in keys]):
            return []

        features = [result['_id']]
        for key in keys:
            features.append(phast_cons[key] if key in phast_cons else np.nan)
        for key in keys:
            features.append(phylop[key] if key in phylop else np.nan)
        # features.extend([phast_cons[key] for key in phast_cons])
        # features.extend([phylop[key] for key in phylop])

        return features


    def _get_uniprot_lines(self, result) -> list:

        if 'snpeff' not in result or 'dbnsfp' not in result:
            return []

        ann = result['snpeff']['ann']
        if isinstance(ann, list):
            hgvs_ps = {ann[ii]['hgvs_p'][2:] for ii in range(len(ann))}
        else:
            hgvs_ps = [ann['hgvs_p'][2:]]
        
        uniprot = result['dbnsfp']['uniprot']
        if isinstance(uniprot, list):
            accs = {uniprot[ii]['acc'].split('-')[0] \
            for ii in range(len(uniprot))}
        else:
            accs = [uniprot['acc'].split('-')[0]]

        ids = []
        for hgvs_p in hgvs_ps:
            for acc in accs:
                a1 = hgvs_p[:3]
                a2 = hgvs_p[-3:]
                loc = hgvs_p[3:-3]

                if a1 not in protein_letters_3to1 \
                    or a2 not in protein_letters_3to1:  # TODO
                    continue

                b1 = protein_letters_3to1[a1]
                b2 = protein_letters_3to1[a2]

                ids.append(f'{acc} {loc} {b1} {b2}')
        
        return ids

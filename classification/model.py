import os.path

import numpy as np
import pandas as pd
from tqdm import tqdm

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.feature_selection import SelectKBest, VarianceThreshold, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, recall_score, precision_score, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight
from validation.cf_matrix import make_confusion_matrix
from xgboost import XGBClassifier

from utils import parser
from validation.cf_matrix import make_confusion_matrix
from imblearn.under_sampling import RandomUnderSampler


config = parser.load_config()
use_polyphen = config['use_polyphen']
classification = config['classification']
cache_dir = config['output']['cache_dir']
classification_dir = config['output']['classification_dir']


def model_trained(SNP_type):
    if SNP_type == 'nsSNP':
        trained_model_fpath = os.path.join(
            classification_dir,
            f'trained_model_{"polyphen" if use_polyphen == "True" else "no_polyphen"}_{classification}_{SNP_type}.model')
    else:
        trained_model_fpath = os.path.join(
            classification_dir,
            f'trained_model_{classification}_{SNP_type}.model')
    return os.path.exists(trained_model_fpath)


class SnpMAPModel:
    
    def __init__(self, dataset_fpath, labels_fpath=None, mode='work', SNP_type='nsSNP'):
        self.data = pd.read_csv(dataset_fpath)
        self.X = self.data.drop(columns='ID').to_numpy()
        self.ids = self.data['ID'].to_list()
        self.model = XGBClassifier()
        self.SNP_type = SNP_type
        if SNP_type == 'nsSNP':
            self.trained_model_fpath = os.path.join(
                classification_dir,
                f'trained_model_{"polyphen" if use_polyphen == "True" else "no_polyphen"}_{classification}_{SNP_type}.model')
        else:
            self.trained_model_fpath = os.path.join(
                classification_dir,
                f'trained_model_{classification}_{SNP_type}.model')
        if mode == 'train':
            self.y = pd.read_csv(labels_fpath).drop(
                columns='ID').to_numpy().flatten()
        if mode == 'work':
            if os.path.exists(self.trained_model_fpath):
                self.model.load_model(self.trained_model_fpath)


    def train(self):

        print('Training the model.')
        
        if self.SNP_type == "sSNP":
            print('undersample')
            rus = RandomUnderSampler(random_state=0, sampling_strategy=0.65)
            self.X, self.y = rus.fit_resample(self.X, self.y)
        scaler = StandardScaler()
        # feature_sel = SelectKBest(f_classif, k='all')
        # variance_sel = VarianceThreshold(threshold=0.00001)
        # imputer = KNNImputer(n_neighbors=32)
        imputer = SimpleImputer(strategy='mean')

        print('Normalizing the features.')
        
        self.X = imputer.fit_transform(self.X)
        # self.X = variance_sel.fit_transform(self.X)
        self.X = scaler.fit_transform(self.X)
        # self.X = feature_sel.fit_transform(self.X, self.y)
        # self.X = scaler.fit_transform(self.X)
        print('Features normalized.')

        label_weights = compute_class_weight(
            'balanced', classes=np.unique(self.y), y=self.y)
        weights = np.asarray([label_weights[j] for j in self.y])

        print(self.X.shape)
        print(self.y.shape)

        objective = 'binary:logistic' \
            if config['classification'] == 'binary' else 'multi:softmax'
        if config['classification'] == 'binary':
            final_model = XGBClassifier(
                n_jobs=-1, objective=objective)
        else:
            final_model = XGBClassifier(
                n_jobs=-1, objective=objective, num_class=4)
        final_model.fit(self.X, self.y, sample_weight=weights)
        final_model.save_model(self.trained_model_fpath)

        train_preds = final_model.predict(self.X)

        if config['classification'] == 'binary':
            print(f'The train AUC is {roc_auc_score(self.y, train_preds)}.')

        make_confusion_matrix(
            confusion_matrix(self.y, train_preds), 
            group_names=[
                'True Negative',
                'False Positive',
                'False Negative',
                'True Positive'],
            categories=['Benign', 'Pathogenic'] \
                if config['classification'] == 'binary' \
                    else ['Benign', 'Likely Benign', 'Likely Pathogenic', 'Pathogenic'],
            cmap='Blues')

        print('Training the model finished.')

    
    def predict(self):
        self.prediction_results = self.model.predict(self.X)
        return {
            'results':{
                self.ids[ii]: self.prediction_results[ii]
                    for ii in range(len(self.prediction_results))},
            'predictions': self.prediction_results}

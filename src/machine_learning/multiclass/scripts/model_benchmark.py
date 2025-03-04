import pandas as pd
import numpy as np
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config

import pickle
np.int = np.int_
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from skopt import BayesSearchCV

script_dir = config.boruta_shap_script
if script_dir not in sys.path:
     sys.path.append(script_dir)
from boruta_shap import UpdatedBorutaShap

import shap
import xgboost
from xgboost import XGBClassifier
from xgboost import plot_importance
import lightgbm
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer, balanced_accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
import seaborn as sn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)
# warnings.filterwarnings("ignore")

import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


seed = 0

dataset = pd.read_csv(config.training_genes, sep="\t")

data = dataset.drop(["Gene"], axis=1)  
data["label_encoded"] = data["label"].map({"most likely": 0, "probable": 1, "least likely": 2})
Y = data["label_encoded"]

data = pd.read_csv(config.cleaned_training_data_ml, header=0, sep=",")
Y = data["label_encoded"]

X = pd.read_csv(config.cleaned_training_data_ml, header=0, sep=",")
X = X.drop(['Gene', "label_encoded"], axis=1)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=seed)

scaler = StandardScaler().fit(X_train)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)

X_train_scaled = X_train_scaled.sparse.to_dense() if hasattr(X_train_scaled, "sparse") else X_train_scaled
X_test_scaled = X_test_scaled.sparse.to_dense() if hasattr(X_test_scaled, "sparse") else X_test_scaled


models = [
    ('XGB', xgboost.XGBClassifier(random_state=seed, objective='multi:softmax', num_class=3, eval_metric='mlogloss')),
    ('LGBM', LGBMClassifier(random_state=seed, verbose=-1)),
    ('CatBoost', CatBoostClassifier(random_seed=seed, verbose=False)),
    ('RF', RandomForestClassifier(random_state=seed)),
    ('LR', LogisticRegression(random_state=seed, multi_class='multinomial')),
    ('SVC', SVC(random_state=seed, probability=True))
]

params = {
    'XGB': config.xgb_parameters,
    'LGBM': config.lgbm_parameters,
    'CatBoost': config.catboost_parameters,
    'RF': config.rf_parameters,
    'LR': config.lr_parameters,
    'SVC': config.svc_parameters
}

inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)

scoring_metrics = {'balanced_accuracy': 'balanced_accuracy', 'f1_macro': 'f1_macro', 'precision_macro': 'precision_macro', 'recall_macro': 'recall_macro'}

model_details = []

for name, model in models:
    if name in ['LR', 'SVC']:
        X_tr, Y_tr = X_train_scaled, Y_train
    else:
        X_tr, Y_tr = X_train, Y_train
    
    search = BayesSearchCV(model, params[name], n_iter=10, cv=inner_cv, n_jobs=-1, random_state=seed)
    search.fit(X_tr, Y_tr)
    
    cv_results = cross_validate(search, X_tr, Y_tr, cv=outer_cv, scoring=scoring_metrics)
    model_details.append((name, search.best_params_, cv_results))

results_df = pd.DataFrame()

for name, best_params, cv_results in model_details:
    for metric in scoring_metrics:
        mean_score = np.mean(cv_results[f'test_{metric}'])
        print(f"{name} - {metric}: {mean_score}")
    print(f"{name} Best Parameters: {best_params}\n")
    df = pd.DataFrame({
        'model': [name],
        'balanced_accuracy': [np.mean(cv_results['test_balanced_accuracy'])],
        'f1_macro': [np.mean(cv_results['test_f1_macro'])],
        'precision_macro': [np.mean(cv_results['test_precision_macro'])],
        'recall_macro': [np.mean(cv_results['test_recall_macro'])]
    })
    results_df = pd.concat([results_df, df], ignore_index=True)

results_df.to_csv(config.ml_eval_metrics, index=False)

best_model_info = max(model_details, key=lambda x: np.mean(x[2]['test_balanced_accuracy']))
best_model_name, best_model_params, _ = best_model_info
print('Top model parameters:')
best_model = [model for name, model in models if name == best_model_name][0]
best_model.set_params(**best_model_params)
print(best_model)

output_dir = config.ml_output_path
best_model.fit(X_train, Y_train)

if best_model_name == 'XGB':
    # Plot feature importance for XGBoost
    plt.figure(figsize=(14, 8))
    xgboost.plot_importance(best_model, importance_type='weight', title='Feature importance (XGBoost - by weight)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/xgboost_feature_importance_weight_no_class_weights.png", bbox_inches='tight', dpi=300)
    plt.close()

elif best_model_name == 'LGBM':
    # Plot feature importance for LightGBM
    plt.figure(figsize=(14, 8))
    lightgbm.plot_importance(best_model, importance_type='split', title='Feature importance (LightGBM - by split)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/lgbm_feature_importance_split_no_class_weights.png", bbox_inches='tight', dpi=300)
    plt.close()

    plt.figure(figsize=(14, 8))
    lightgbm.plot_importance(best_model, importance_type='gain', title='Feature importance (LightGBM - by gain)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/lgbm_feature_importance_gain_no_class_weights.png", bbox_inches='tight', dpi=300)
    plt.close()

elif best_model_name == 'CatBoost':
    # Plot feature importance for CatBoost
    feature_importances = best_model.get_feature_importance(type="FeatureImportance")
    feature_names = X_train.columns
    plt.figure(figsize=(14, 8))
    plt.bar(range(len(feature_importances)), feature_importances)
    plt.xticks(ticks=range(len(feature_names)), labels=feature_names, rotation='vertical')
    plt.title('Feature importance (CatBoost)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/catboost_feature_importance_no_class_weights.png", bbox_inches='tight', dpi=300)
    plt.close()

else:
    print(f"No feature importance plotting is available for the model {best_model_name}.")



#Feature_Selector = UpdatedBorutaShap(importance_measure="shap", classification=False)
print('Boruta Shap Feature Importance:')
Feature_Selector = UpdatedBorutaShap(importance_measure="shap", classification=False)
Feature_Selector.fit(X=X_train, y=Y_train, n_trials=100, random_state=seed)

Feature_Selector.plot(which_features="all", X_rotation=90, figsize=(18,8))
fig = plt.gcf()  
fig.savefig(config.boruta_shap_plot, dpi=300, bbox_inches='tight')  
plt.close(fig) 

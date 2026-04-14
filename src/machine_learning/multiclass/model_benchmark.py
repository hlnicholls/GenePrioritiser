import pandas as pd
import numpy as np
import sys
import os
import resource
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
CONFIG_DIR = os.path.join(ROOT_DIR, 'config')
if CONFIG_DIR not in sys.path:
    sys.path.insert(0, CONFIG_DIR)
import config

# ── Resource limits (set in config/config.py) ───────────────────────────────────────
MAX_CORES = config.MAX_CORES
MAX_RAM_GB = config.MAX_RAM_GB
MAX_RAM_BYTES = MAX_RAM_GB * 1024 ** 3

try:
    resource.setrlimit(resource.RLIMIT_AS, (MAX_RAM_BYTES, MAX_RAM_BYTES))
    print(f"[config] Memory limit set to {MAX_RAM_GB} GB")
except Exception as e:
    print(f"[config] Warning: could not set memory limit: {e}")

os.environ['OMP_NUM_THREADS'] = str(MAX_CORES)
os.environ['MKL_NUM_THREADS'] = str(MAX_CORES)
os.environ['OPENBLAS_NUM_THREADS'] = str(MAX_CORES)
print(f"[config] CPU threads limited to {MAX_CORES} cores")

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

print("\n[data] Loading training data...")
dataset = pd.read_csv(config.training_genes, sep="\t")

data = dataset.drop(["Gene"], axis=1)  
data["label_encoded"] = data["label"].map({"most likely": 0, "probable": 1, "least likely": 2})
Y = data["label_encoded"]

data = pd.read_csv(config.cleaned_training_data_ml, header=0, sep=",")
Y = data["label_encoded"]

X = pd.read_csv(config.cleaned_training_data_ml, header=0, sep=",")
X = X.drop(['Gene', "label_encoded"], axis=1)
print(f"[data] Loaded {X.shape[0]} samples, {X.shape[1]} features")
print(f"[data] Label distribution: {Y.value_counts().to_dict()}")
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=seed)
print(f"[data] Train: {X_train.shape[0]} samples | Test: {X_test.shape[0]} samples")

scaler = StandardScaler().fit(X_train)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)

X_train_scaled = X_train_scaled.sparse.to_dense() if hasattr(X_train_scaled, "sparse") else X_train_scaled
X_test_scaled = X_test_scaled.sparse.to_dense() if hasattr(X_test_scaled, "sparse") else X_test_scaled


models = [
    ('XGB', xgboost.XGBClassifier(random_state=seed, objective='multi:softmax', num_class=3, eval_metric='mlogloss', nthread=1)),
    ('LGBM', LGBMClassifier(random_state=seed, verbose=-1, n_jobs=1)),
    ('CatBoost', CatBoostClassifier(random_seed=seed, verbose=False, thread_count=1)),
    ('RF', RandomForestClassifier(random_state=seed, n_jobs=1)),
    ('LR', LogisticRegression(random_state=seed, multi_class='multinomial', n_jobs=1)),
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
print(f"\n[benchmark] Starting benchmark for {len(models)} models...")

for i, (name, model) in enumerate(models, 1):
    if name in ['LR', 'SVC']:
        X_tr, Y_tr = X_train_scaled, Y_train
    else:
        X_tr, Y_tr = X_train, Y_train

    print(f"\n[{i}/{len(models)}] {name}: running Bayesian hyperparameter search (n_iter=10, cv=5)...")
    search = BayesSearchCV(model, params[name], n_iter=10, cv=inner_cv, n_jobs=MAX_CORES, random_state=seed)
    search.fit(X_tr, Y_tr)
    print(f"  Best params: {dict(search.best_params_)}")

    print(f"  Running outer cross-validation (cv=5)...")
    cv_results = cross_validate(search, X_tr, Y_tr, cv=outer_cv, scoring=scoring_metrics, n_jobs=MAX_CORES)
    mean_ba = np.mean(cv_results['test_balanced_accuracy'])
    mean_f1 = np.mean(cv_results['test_f1_macro'])
    print(f"  CV balanced_accuracy: {mean_ba:.4f} | f1_macro: {mean_f1:.4f}")
    model_details.append((name, search.best_params_, cv_results))

print("\n[results] Summarising benchmark results...")
results_df = pd.DataFrame()

for name, best_params, cv_results in model_details:
    for metric in scoring_metrics:
        mean_score = np.mean(cv_results[f'test_{metric}'])
        print(f"  {name} - {metric}: {mean_score:.4f}")
    print(f"  {name} Best Parameters: {best_params}\n")
    df = pd.DataFrame({
        'model': [name],
        'balanced_accuracy': [np.mean(cv_results['test_balanced_accuracy'])],
        'f1_macro': [np.mean(cv_results['test_f1_macro'])],
        'precision_macro': [np.mean(cv_results['test_precision_macro'])],
        'recall_macro': [np.mean(cv_results['test_recall_macro'])]
    })
    results_df = pd.concat([results_df, df], ignore_index=True)

os.makedirs(os.path.dirname(config.ml_eval_metrics), exist_ok=True)
results_df.to_csv(config.ml_eval_metrics, index=False)
print(f"[results] Metrics saved to {config.ml_eval_metrics}")

best_model_info = max(model_details, key=lambda x: np.mean(x[2]['test_balanced_accuracy']))
best_model_name, best_model_params, _ = best_model_info
print(f"\n[best model] Top model: {best_model_name}")
print(f"[best model] Best parameters: {dict(best_model_params)}")
best_model = [model for name, model in models if name == best_model_name][0]
best_model.set_params(**best_model_params)
print(best_model)

output_dir = config.ml_output_path
os.makedirs(output_dir, exist_ok=True)
print(f"[best model] Fitting {best_model_name} on full training set...")
best_model.fit(X_train, Y_train)
print(f"[best model] Fit complete.")

print(f"\n[feature importance] Plotting feature importance for {best_model_name}...")
if best_model_name == 'XGB':
    plt.figure(figsize=(14, 8))
    xgboost.plot_importance(best_model, importance_type='weight', title='Feature importance (XGBoost - by weight)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/xgboost_feature_importance_weight_no_class_weights.png", bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  Saved XGBoost feature importance plot.")

elif best_model_name == 'LGBM':
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
    print(f"  Saved LightGBM feature importance plots.")

elif best_model_name == 'CatBoost':
    feature_importances = best_model.get_feature_importance(type="FeatureImportance")
    feature_names = X_train.columns
    plt.figure(figsize=(14, 8))
    plt.bar(range(len(feature_importances)), feature_importances)
    plt.xticks(ticks=range(len(feature_names)), labels=feature_names, rotation='vertical')
    plt.title('Feature importance (CatBoost)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/catboost_feature_importance_no_class_weights.png", bbox_inches='tight', dpi=300)
    plt.close()
    print(f"  Saved CatBoost feature importance plot.")

else:
    print(f"  No feature importance plotting available for {best_model_name}.")



print("\n[boruta] Running BorutaSHAP feature selection (n_trials=100) — this may take a while...")
Feature_Selector = UpdatedBorutaShap(importance_measure="shap", classification=False)
Feature_Selector.fit(X=X_train, y=Y_train, n_trials=100, random_state=seed)
print("[boruta] BorutaSHAP complete.")

print("[boruta] Saving BorutaSHAP plot...")
os.makedirs(os.path.dirname(config.boruta_shap_plot), exist_ok=True)
Feature_Selector.plot(which_features="all", X_rotation=90, figsize=(18,8))
fig = plt.gcf()  
fig.savefig(config.boruta_shap_plot, dpi=300, bbox_inches='tight')  
plt.close(fig)
print(f"[boruta] Plot saved to {config.boruta_shap_plot}")
print("\n[done] model_benchmark.py complete.") 

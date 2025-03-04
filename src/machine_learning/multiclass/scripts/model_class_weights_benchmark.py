import pandas as pd
import numpy as np
import pickle
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, Input
from scikeras.wrappers import KerasClassifier
import sys
script_dir = config.eda_script_path
sys.path.append(script_dir)
from updated_MissForest_Algorithm import MissForest
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from skopt import BayesSearchCV
import xgboost
from xgboost import XGBClassifier, plot_importance
import lightgbm
from lightgbm import LGBMClassifier
from catboost import CatBoostClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import make_scorer, balanced_accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
from sklearn.utils.class_weight import compute_class_weight
script_dir = config.boruta_shap_script
if script_dir not in sys.path:
     sys.path.append(script_dir)
from boruta_shap import UpdatedBorutaShap

import shap
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore")

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

class_weights = compute_class_weight('balanced', classes=np.unique(Y_train), y=Y_train)
weight_dict = dict(zip(np.unique(Y_train), class_weights))

def create_nn_model():
    model = Sequential()
    model.add(Input(shape=(X_train.shape[1],)))  # Explicit Input layer
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(3, activation='softmax'))
    model.compile(optimizer='adam', loss='sparse_categorical_crossentropy', metrics=['accuracy'])
    return model

nn_model = KerasClassifier(model=create_nn_model, epochs=50, batch_size=32, verbose=0)

models = [
    ('XGB', XGBClassifier(random_state=seed, objective='multi:softmax', num_class=3, eval_metric='mlogloss')),
    ('LGBM', LGBMClassifier(random_state=seed, verbose=-1, class_weight='balanced')),
    ('CatBoost', CatBoostClassifier(random_seed=seed, verbose=False, auto_class_weights='Balanced')),
    ('RF', RandomForestClassifier(random_state=seed, class_weight='balanced')),
    ('LR', LogisticRegression(random_state=seed, multi_class='multinomial', class_weight='balanced')),
    ('SVC', SVC(random_state=seed, probability=True, class_weight='balanced')),
    ('NN', nn_model)
]

params = {
    'XGB': config.xgb_parameters,
    'LGBM': config.lgbm_parameters,
    'CatBoost': config.catboost_parameters,
    'RF': config.rf_parameters,
    'LR': config.lr_parameters,
    'SVC': config.svc_parameters,
    'NN': config.nn_parameters
}


inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)

scoring_metrics = {'balanced_accuracy': 'balanced_accuracy', 'f1_macro': 'f1_macro', 'precision_macro': 'precision_macro', 'recall_macro': 'recall_macro'}

model_details = []
fold_balanced_accuracies = []

for name, model in models:
    if name in ['LR', 'SVC', 'NN']:
        X_tr, Y_tr = X_train_scaled, Y_train
    else:
        X_tr, Y_tr = X_train, Y_train
    
    # Assuming `search` is the BayesSearchCV object:
    search = BayesSearchCV(model, params[name], n_iter=10, cv=inner_cv, n_jobs=-1, random_state=seed)
    search.fit(X_tr, Y_tr)
    
    cv_results = cross_validate(search, X_tr, Y_tr, cv=outer_cv, scoring=scoring_metrics)
    model_details.append((name, search.best_params_, cv_results))

    fold_balanced_accuracies.append(cv_results['test_balanced_accuracy'])

# Create DataFrame to store balanced accuracy for each fold for each model
fold_balanced_accuracies_df = pd.DataFrame(fold_balanced_accuracies, index=[name for name, _ in models])
fold_balanced_accuracies_df.columns = [f'Fold_{i+1}' for i in range(fold_balanced_accuracies_df.shape[1])]

# Save the DataFrame as a CSV file
fold_balanced_accuracies_df.to_csv(config.balanced_accuracy_per_fold, index=True)

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

results_df.to_csv(config.ml_eval_metrics_class_weighted, index=False)

# Find the top 3 performing models based on balanced accuracy
top_3_models = sorted(model_details, key=lambda x: np.mean(x[2]['test_balanced_accuracy']), reverse=True)[:3]

# Create the VotingClassifier using the top 3 models
# voting_clf = VotingClassifier(estimators=[(name, models[[m[0] for m in models].index(name)][1]) for name, _, _ in top_3_models], voting='soft')
# voting_clf.fit(X_train, Y_train)

# # Save the VotingClassifier model
# filename = config.best_voting_model_pkl_file
# with open(filename, 'wb') as file:
#     pickle.dump(voting_clf, file)

# Fit and save the best individual model
best_model_info = top_3_models[0]
best_model_name, best_model_params, _ = best_model_info
best_model = models[[m[0] for m in models].index(best_model_name)][1]
best_model.set_params(**best_model_params)
best_model.fit(X_train, Y_train)

output_dir = config.ml_output_path

if best_model_name == 'XGB':
    # Plot feature importance for XGBoost
    plt.figure(figsize=(14, 8))
    xgboost.plot_importance(best_model, importance_type='weight', title='Feature importance (XGBoost - by weight)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/xgboost_feature_importance_weight.png", bbox_inches='tight', dpi=300)
    plt.close()

elif best_model_name == 'LGBM':
    # Plot feature importance for LightGBM
    plt.figure(figsize=(14, 8))
    lightgbm.plot_importance(best_model, importance_type='split', title='Feature importance (LightGBM - by split)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/lgbm_feature_importance_split_class_weighted.png", bbox_inches='tight', dpi=300)
    plt.close()

    plt.figure(figsize=(14, 8))
    lightgbm.plot_importance(best_model, importance_type='gain', title='Feature importance (LightGBM - by gain)')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/lgbm_feature_importance_gaint_class_weighted.png", bbox_inches='tight', dpi=300)
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
    plt.savefig(f"{output_dir}/catboost_feature_importancet_class_weighted.png", bbox_inches='tight', dpi=300)
    plt.close()

else:
    print(f"No feature importance plotting is available for the model {best_model_name}.")

Feature_Selector = UpdatedBorutaShap(importance_measure="shap", classification=False)
Feature_Selector.fit(X=X_train, y=Y_train, n_trials=100, random_state=seed)
fig = plt.gcf()  
fig.savefig(config.boruta_shap_plot_class_weighted, dpi=300, bbox_inches='tight')  
plt.close(fig) 


print('Saving best model...')
best_model.fit(X, Y)
filename = config.best_model_pkl_file
with open(filename, 'wb') as file:
    pickle.dump(best_model, file)

# SHAP analysis for the best model
explainer = shap.TreeExplainer(best_model)
shap_values = explainer.shap_values(X_test)

for class_idx in range(len(shap_values)):
    plt.figure()
    shap.summary_plot(shap_values[class_idx], X_test, show=False)
    plt.grid(True)
    plt.title(f"SHAP Summary Plot for Class {class_idx}")
    plt.savefig(f"{output_dir}/shap_summary_class_{class_idx}.png", bbox_inches='tight', dpi=300)
    plt.close()

# Generate predictions and plot confusion matrix
Y_pred = best_model.predict(X_test)
conf_matrix = confusion_matrix(Y_test, Y_pred)

labels = ["Most Likely", "Probable", "Least Likely"]

plt.figure(figsize=(10, 8))
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=labels, yticklabels=labels, annot_kws={"size": 16})
plt.xlabel('Predicted Label')
plt.ylabel('True Label')
plt.title(f'Confusion Matrix - {best_model_name}')
plt.tight_layout()
plt.savefig(f"{output_dir}/confusion_matrix_{best_model_name}.png", bbox_inches='tight', dpi=300)
plt.close()

print("Model training and evaluation completed. Outputs have been saved.")

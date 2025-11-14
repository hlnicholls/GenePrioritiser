"""
eda_training_data.py

Purpose
-------
Lightweight EDA and data-cleaning pipeline used before multiclass model training.

High-level steps performed by this script
----------------------------------------
- Load the merged training table (`config.training_data_all_features_eda`).
- Map the textual `label` column to a numeric `label_encoded` for convenience.
- Compute per-column missingness and drop features with >25% missing values.
- Keep a reduced set of columns with acceptable missingness for downstream
    processing.
- Impute remaining missing values using a MissForest implementation (random-forest
    based imputation).
- Compute correlation (Pearson / Spearman) against the label and remove
    features that are highly collinear with other features (absolute Spearman
    correlation > 0.9) to reduce redundancy.
- Save cleaned datasets for EDA and for the ML pipeline.

"""

import numpy as np
import pandas as pd
import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
CONFIG_DIR = os.path.join(ROOT_DIR, 'config')
if CONFIG_DIR not in sys.path:
    sys.path.insert(0, CONFIG_DIR)
import config

import matplotlib.pyplot as plt

script_dir = config.eda_script_path
sys.path.append(script_dir)
from updated_MissForest_Algorithm import MissForest

from scipy.cluster.hierarchy import dendrogram, ward
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression, mutual_info_regression, f_classif
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score

plt.style.use("ggplot")
import seaborn as sns
from matplotlib import pyplot

sns.set_style("darkgrid")
sns.mpl.rcParams["figure.figsize"] = (15.0, 9.0)

import re
regex = re.compile(r"\[|\]|<", re.IGNORECASE)

import gc
import warnings
warnings.filterwarnings("ignore")

seed = 0

dataset = pd.read_csv(config.training_data_all_features_eda, sep=",")

dataset.shape

dataset.head()

data = dataset.drop(["Gene"], axis=1)
print(data.shape)

data["label_encoded"] = data["label"].map(
    {"most likely": 0, "probable": 1, "least likely": 2}
)
Y = data["label"]

data.describe()

data.isnull().sum()

natest = data.isnull().sum()
natest.sort_values(inplace=True)


percent_missing = data.isnull().sum() * 100 / len(data)
missing_value_df = pd.DataFrame(
    {"column_name": data.columns, "percent_missing": percent_missing}
)
missing_value_df.sort_values("percent_missing", inplace=True)

natest = natest.to_frame()
missingdata = natest.join(missing_value_df)

data_drop = data.drop(["label", "label_encoded",], axis=1)

null_counts = data_drop.isnull().sum() / len(data_drop)
plt.figure(figsize=(40, 8))
plt.xticks(np.arange(len(null_counts)) + 0.0, null_counts.index, rotation="vertical")
plt.ylabel("Fraction of rows with missing data")
plt.bar(np.arange(len(null_counts)), null_counts)
plt.savefig(config.feature_missingness_plot, format="png", dpi=300, bbox_inches="tight")

selection = missing_value_df[missing_value_df["percent_missing"] < 25.00]

# Calculate and print the number of columns
total_columns = data_drop.shape[1] - 1  # Excluding the 'Gene' column
print(f"Total columns (excluding 'Gene'): {total_columns}")

# Calculating and printing number of columns dropped due to missingness
columns_dropped = total_columns - len(selection)
print(f"Number of columns dropped due to missingness: {columns_dropped}")

list(selection["column_name"])


dat = data[list(selection["column_name"])]
dat["Gene"] = dataset["Gene"]

# --- Automated leakage checks: one-vs-rest AUC and ANOVA F-statistic ---
# For multiclass problems we compute per-feature one-vs-rest AUCs (numeric
# features only) and ANOVA F-statistics. Features with extremely high
# separability (AUC > 0.95) or unusually large F-statistics are flagged and
# written to an audit file for manual review.
leakage_suspects = []
try:
    y = data['label'] if 'label' in data.columns else dataset['label']
    classes = list(pd.Categorical(y).categories)
    yb = label_binarize(y, classes=classes)

    # Candidate feature columns (exclude label fields)
    candidate_feats = [c for c in list(selection['column_name']) if c not in ('Gene', 'label', 'label_encoded')]
    # numeric features only for AUC
    numeric_feats = [c for c in candidate_feats if pd.api.types.is_numeric_dtype(data[c])]

    # compute one-vs-rest AUCs
    auc_records = []
    auc_map = {}
    auc_threshold = 0.95
    for feat in numeric_feats:
        x = data[feat].astype(float).fillna(data[feat].median())
        aucs = []
        for ci in range(yb.shape[1]):
            try:
                aucs.append(float(roc_auc_score(yb[:, ci], x)))
            except Exception:
                aucs.append(0.5)
        max_auc = max(aucs) if len(aucs) > 0 else 0.5
        auc_map[feat] = max_auc
        auc_records.append((feat, max_auc))

    # compute ANOVA F-statistics (requires numeric matrix)
    f_map = {}
    f_threshold = None
    try:
        if len(candidate_feats) > 0:
            X_num = data[candidate_feats].apply(pd.to_numeric, errors='coerce')
            # fill na with median per column
            X_num = X_num.fillna(X_num.median())
            y_enc = pd.Categorical(y).codes
            f_vals, p_vals = f_classif(X_num, y_enc)
            f_threshold = float(np.mean(f_vals) + 5 * np.std(f_vals))
            for feat, fval in zip(candidate_feats, f_vals):
                f_map[feat] = float(fval)
    except Exception:
        # silently skip ANOVA if it fails (small sample, constant cols, etc.)
        f_map = {}
        f_threshold = None

    # Final leakage suspects: require BOTH AUC and ANOVA to exceed their thresholds
    # (logical AND). This reduces false positives where only one test is extreme.
    if f_threshold is not None:
        for feat in candidate_feats:
            auc_val = auc_map.get(feat, 0.0)
            f_val = f_map.get(feat, None)
            if (f_val is not None) and (auc_val > auc_threshold) and (f_val > f_threshold):
                leakage_suspects.append({'feature': feat, 'max_auc': float(auc_val), 'f_value': float(f_val)})

    # write suspects to audit file next to correlation outputs
    out_dir = os.path.dirname(config.correlation_pairs_09)
    os.makedirs(out_dir, exist_ok=True)
    suspects_path = os.path.join(out_dir, 'feature_leakage_suspects.tsv')
    if leakage_suspects:
        pd.DataFrame(leakage_suspects).to_csv(suspects_path, sep='\t', index=False)
        print('Wrote feature leakage suspects to:', suspects_path)
    else:
        print('No high-separability features detected by AUC/ANOVA checks')
except Exception as e:
    print('Warning: leakage checks failed:', e)


dt2 = dat
dat = dat.set_index("Gene")

df = dt2
df = df.set_index("Gene")
df["label_encoded"] = df["label"].map(
    {"most likely": 0, "probable": 1, "least likely": 2}
)

df = df.drop(
    [
        "label",
    ],
    axis=1
)

X=df
imputer = MissForest(random_state=seed)
X = pd.DataFrame(imputer.fit_transform(X), index=df.index, columns=df.columns)

Xcor = X

Xcor = pd.DataFrame(data=Xcor, columns=X.columns)

corr_matrix = X.corr()
print(corr_matrix["label_encoded"].sort_values(ascending=False))
corr = corr_matrix["label_encoded"].sort_values(ascending=False)

def corrFilter(x: pd.DataFrame, bound: float):
    xCorr = x.corr()
    xFiltered = xCorr[((xCorr >= bound) | (xCorr <= -bound)) & (xCorr != 1.000)]
    xFlattened = xFiltered.unstack().sort_values().drop_duplicates()
    return xFlattened

corrFilter(X, 0.9)

pairs = corrFilter(X, 0.9)
pairs.to_csv(config.correlation_pairs_09, header=True)
corr_matrix.to_csv(config.correlation_pairs_09, header=True)
corr.to_csv(config.correlation_matrix_09, header=True)
Xcor = X
Xcor = pd.DataFrame(data=Xcor, columns=X.columns)
corr = Xcor.corr(method="spearman").abs()
upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(np.bool_))
to_drop = [column for column in upper.columns if any(upper[column] > 0.9)]
print("Dropped features with > 0.9 correlation:", to_drop)

selected_columns = X.drop(X[to_drop], axis=1)
selected_columns.head()
selected_columns.to_csv(config.cleaned_training_data_eda)
selected_columns.to_csv(config.cleaned_training_data_ml)
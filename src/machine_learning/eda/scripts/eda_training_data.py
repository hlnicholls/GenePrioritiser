import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config

import matplotlib.pyplot as plt

script_dir = config.eda_script_path
sys.path.append(script_dir)
from updated_MissForest_Algorithm import MissForest

from scipy.cluster.hierarchy import dendrogram, ward
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression, mutual_info_regression
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler

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
#dataset = dataset.rename({"IPA_BP": "IPA_BP_annotation"}, axis=1)
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

#missingdata.to_csv("traininggenes_features_with_NA.csv")

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

#X = MinMaxScaler().fit_transform(df)
X=df
imputer = MissForest(random_state=seed)
X = pd.DataFrame(imputer.fit_transform(X), index=df.index, columns=df.columns)

Xcor = X

Xcor = pd.DataFrame(data=Xcor, columns=X.columns)

# Heatmap:
# corr = Xcor.corr()
# f, ax = plt.subplots(figsize=(60, 50))

# cmap = sns.diverging_palette(10, 275, as_cmap=True)
# sns.set(font_scale=3)

# htmap = sns.heatmap(
#     corr,
#     cmap=cmap,
#     square=True,
#     xticklabels=True,
#     yticklabels=True,
#     linewidths=0.5,
#     ax=ax,
# )

# figure = htmap.get_figure()
# figure.savefig("/data/WHRI-Bioinformatics/Nicholls/BP_Gene_Prioritisation/BP_Gene_Prior_v2/2_machine_learning/EDA/Output/correlation_heatmap.png", format="png", dpi=300, bbox_inches="tight")

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
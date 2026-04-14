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
from scipy.stats import mannwhitneyu as _mannwhitneyu
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from sklearn.feature_selection import f_regression, mutual_info_regression, f_classif
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_auc_score
from statsmodels.stats.outliers_influence import variance_inflation_factor

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


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


EDA_OUT_DIR = os.path.dirname(config.correlation_pairs_09)
os.makedirs(EDA_OUT_DIR, exist_ok=True)


def eda_out(name):
    return os.path.join(EDA_OUT_DIR, name)

dataset = pd.read_csv(config.training_data_all_features_eda, sep=",")

dataset.shape

dataset.head()

data = dataset.drop(["Gene"], axis=1)
print(data.shape)

data["label_encoded"] = data["label"].map(
    {"most likely": 0, "probable": 1, "least likely": 2}
)
Y = data["label"]

# Publication plot 1: class balance
class_counts = Y.value_counts().reindex(["most likely", "probable", "least likely"], fill_value=0)
plt.figure(figsize=(8, 5))
sns.barplot(x=class_counts.index, y=class_counts.values, palette="Set2")
plt.title("Training Class Balance")
plt.xlabel("Class")
plt.ylabel("Gene count")
plt.xticks(rotation=15)
plt.tight_layout()
class_balance_plot = eda_out("class_balance.png")
plt.savefig(class_balance_plot, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved class balance plot to: {class_balance_plot}")

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
ensure_parent_dir(config.feature_missingness_plot)
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
auc_map = {}
f_map = {}
f_threshold = None
auc_threshold = 0.95
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

# Publication plot 2: missingness before vs after imputation/filtering
missing_before = data_drop.isnull().mean() * 100
missing_after = X.drop(columns=["label_encoded"], errors="ignore").isnull().mean() * 100
missing_compare = pd.DataFrame(
    {
        "feature": missing_before.index,
        "missing_before_pct": missing_before.values,
        "missing_after_pct": missing_after.reindex(missing_before.index, fill_value=0.0).values,
    }
)
missing_compare_path = eda_out("missingness_before_after.tsv")
missing_compare.to_csv(missing_compare_path, sep="\t", index=False)

plt.figure(figsize=(9, 5))
plt.hist(missing_before.values, bins=20, alpha=0.7, label="Before imputation")
plt.hist(missing_after.values, bins=20, alpha=0.7, label="After imputation")
plt.xlabel("Missingness (%)")
plt.ylabel("Number of features")
plt.title("Feature Missingness Before vs After Imputation")
plt.legend()
plt.tight_layout()
missing_compare_plot = eda_out("missingness_before_after.png")
plt.savefig(missing_compare_plot, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved missingness comparison table to: {missing_compare_path}")
print(f"Saved missingness comparison plot to: {missing_compare_plot}")

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

corrFilter(X, 0.98)

pairs = corrFilter(X, 0.98)
ensure_parent_dir(config.correlation_pairs_09)
pairs.to_csv(config.correlation_pairs_09, header=True)
corr_matrix.to_csv(config.correlation_matrix_09, header=True)
print(f"Saved correlation pairs to: {config.correlation_pairs_09}")
print(f"Saved full correlation matrix to: {config.correlation_matrix_09}")
Xcor = X
Xcor = pd.DataFrame(data=Xcor, columns=X.columns)
corr = Xcor.corr(method="spearman").abs()
upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(np.bool_))
corr_threshold = 0.98
to_drop = [column for column in upper.columns if any(upper[column] > corr_threshold)]
print("Dropped features with > 0.98 correlation:", to_drop)

# Save explicit high-correlation feature pairs for interpretability.
pair_table = (
    upper.stack()
    .reset_index()
    .rename(columns={"level_0": "feature_a", "level_1": "feature_b", 0: "abs_spearman_corr"})
)
pair_table = pair_table[pair_table["abs_spearman_corr"] >= corr_threshold]
pair_table = pair_table[(pair_table["feature_a"] != "label_encoded") & (pair_table["feature_b"] != "label_encoded")]
pair_table = pair_table.sort_values("abs_spearman_corr", ascending=False)
pair_table.to_csv(config.correlation_pairs_09, index=False)
print(f"Saved high-correlation pair table (|rho| >= {corr_threshold}) to: {config.correlation_pairs_09}")

# For each dropped feature, report which feature it is most correlated with.
retained_features = [c for c in corr.columns if c not in to_drop and c != "label_encoded"]
dropped_corr_rows = []
for dropped_feature in to_drop:
    if dropped_feature == "label_encoded" or dropped_feature not in corr.columns:
        continue

    # strongest correlation among retained features
    best_partner = None
    best_corr = np.nan
    if retained_features:
        s = corr.loc[dropped_feature, retained_features]
        if len(s) > 0:
            best_partner = s.idxmax()
            best_corr = float(s.max())

    # fallback: strongest correlation among all other non-label features
    if best_partner is None:
        all_candidates = [c for c in corr.columns if c not in {dropped_feature, "label_encoded"}]
        if all_candidates:
            s2 = corr.loc[dropped_feature, all_candidates]
            best_partner = s2.idxmax()
            best_corr = float(s2.max())

    dropped_corr_rows.append(
        {
            "dropped_feature": dropped_feature,
            "strongest_partner_feature": best_partner,
            "abs_spearman_corr": best_corr,
            "partner_retained": bool(best_partner in retained_features) if best_partner is not None else False,
        }
    )

out_dir = os.path.dirname(config.correlation_pairs_09)
os.makedirs(out_dir, exist_ok=True)
dropped_map_path = os.path.join(out_dir, "dropped_feature_correlation_map.tsv")
pd.DataFrame(dropped_corr_rows).to_csv(dropped_map_path, sep="\t", index=False)
print(f"Saved dropped-feature correlation map to: {dropped_map_path}")

# Plot heatmap of dropped features vs their strongest correlated partner features.
partner_features = [
    r["strongest_partner_feature"]
    for r in dropped_corr_rows
    if r.get("strongest_partner_feature") is not None
]
partner_features = list(dict.fromkeys(partner_features))
if len(to_drop) > 0 and len(partner_features) > 0:
    heatmap_mat = corr.loc[to_drop, partner_features]
    plt.figure(figsize=(max(10, len(partner_features) * 0.45), max(8, len(to_drop) * 0.28)))
    sns.heatmap(heatmap_mat, cmap="viridis", vmin=0.0, vmax=1.0)
    plt.title(f"Abs Spearman Correlation: Dropped Features vs Strongest Partner (>|{corr_threshold}|)")
    plt.xlabel("Strongest correlated feature")
    plt.ylabel("Dropped feature")
    plt.tight_layout()
    corr_heatmap_path = os.path.join(out_dir, "dropped_feature_correlation_heatmap.png")
    plt.savefig(corr_heatmap_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved dropped-feature correlation heatmap to: {corr_heatmap_path}")

selected_columns = X.drop(X[to_drop], axis=1)
selected_columns.head()

# Publication plot 3: leakage diagnostic scatter (AUC vs ANOVA F)
if len(auc_map) > 0 and len(f_map) > 0:
    leakage_df = pd.DataFrame(
        {
            "feature": list(set(list(auc_map.keys()) + list(f_map.keys()))),
        }
    )
    leakage_df["max_auc_ovr"] = leakage_df["feature"].map(auc_map)
    leakage_df["anova_f"] = leakage_df["feature"].map(f_map)
    leakage_df = leakage_df.dropna(subset=["max_auc_ovr", "anova_f"])
    leakage_df["is_flagged"] = (
        (leakage_df["max_auc_ovr"] > auc_threshold)
        & (leakage_df["anova_f"] > (f_threshold if f_threshold is not None else np.inf))
    )

    leakage_scatter_path = eda_out("feature_leakage_auc_vs_f.png")
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=leakage_df,
        x="max_auc_ovr",
        y="anova_f",
        hue="is_flagged",
        palette={False: "#4C78A8", True: "#E45756"},
    )
    if f_threshold is not None:
        plt.axhline(f_threshold, color="black", linestyle="--", linewidth=1, label="F threshold")
    plt.axvline(auc_threshold, color="gray", linestyle="--", linewidth=1, label="AUC threshold")
    plt.title("Leakage Diagnostics: One-vs-Rest AUC vs ANOVA F")
    plt.xlabel("Max one-vs-rest AUC")
    plt.ylabel("ANOVA F-statistic")
    plt.tight_layout()
    plt.savefig(leakage_scatter_path, dpi=300, bbox_inches="tight")
    plt.close()
    leakage_df.to_csv(eda_out("feature_leakage_auc_vs_f.tsv"), sep="\t", index=False)
    print(f"Saved leakage diagnostic scatter to: {leakage_scatter_path}")

# Publication plot 4: PCA projection by class (on retained features)
X_retained = selected_columns.drop(columns=["label_encoded"], errors="ignore")
y_codes = selected_columns["label_encoded"].astype(int)
label_map = {0: "most likely", 1: "probable", 2: "least likely"}
y_labels = y_codes.map(label_map)

if X_retained.shape[1] >= 2:
    pca = PCA(n_components=2, random_state=seed)
    X_pca = pca.fit_transform(X_retained)
    pca_df = pd.DataFrame(
        {
            "PC1": X_pca[:, 0],
            "PC2": X_pca[:, 1],
            "label": y_labels.values,
        },
        index=X_retained.index,
    )
    pca_df.to_csv(eda_out("pca_projection.tsv"), sep="\t", index=True)
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="label", alpha=0.8)
    plt.title(
        f"PCA of Retained Features (PC1 {pca.explained_variance_ratio_[0]*100:.1f}%, "
        f"PC2 {pca.explained_variance_ratio_[1]*100:.1f}%)"
    )
    plt.tight_layout()
    pca_plot_path = eda_out("pca_by_class.png")
    plt.savefig(pca_plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved PCA plot to: {pca_plot_path}")

# Publication plot 5: clustermap of retained feature correlations
if X_retained.shape[1] >= 2:
    corr_retained = X_retained.corr(method="spearman")
    corr_retained.to_csv(eda_out("retained_feature_spearman_matrix.tsv"), sep="\t")
    cg = sns.clustermap(
        corr_retained,
        cmap="vlag",
        center=0,
        figsize=(12, 12),
        xticklabels=False,
        yticklabels=False,
    )
    cg.fig.suptitle("Clustermap of Retained Feature Correlations (Spearman)", y=1.02)
    clustermap_path = eda_out("retained_feature_clustermap.png")
    cg.savefig(clustermap_path, dpi=300, bbox_inches="tight")
    plt.close(cg.fig)
    print(f"Saved retained-feature clustermap to: {clustermap_path}")

# Publication plot 6: univariate effect-size panel (ANOVA F on retained features)
retained_f_rows = []
for feat in X_retained.columns:
    if feat in f_map:
        retained_f_rows.append({"feature": feat, "anova_f": f_map[feat]})
retained_f_df = pd.DataFrame(retained_f_rows)
if not retained_f_df.empty:
    retained_f_df = retained_f_df.sort_values("anova_f", ascending=False)
retained_f_df.to_csv(eda_out("retained_feature_anova_f.tsv"), sep="\t", index=False)

if not retained_f_df.empty:
    top_n = min(20, len(retained_f_df))
    top_f = retained_f_df.head(top_n).iloc[::-1]
    plt.figure(figsize=(10, max(6, top_n * 0.35)))
    plt.barh(top_f["feature"], top_f["anova_f"])
    plt.xlabel("ANOVA F-statistic")
    plt.ylabel("Feature")
    plt.title("Top Retained Features by Univariate Class-Separation (ANOVA F)")
    plt.tight_layout()
    effect_plot_path = eda_out("retained_feature_effect_size_anovaF.png")
    plt.savefig(effect_plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved ANOVA effect-size plot to: {effect_plot_path}")

# Publication plot 7: feature distributions by class (features ranked by ANOVA F class-separation)
if not retained_f_df.empty:
    top_features = retained_f_df.head(min(12, len(retained_f_df)))["feature"].tolist()
    class_order = ["most likely", "probable", "least likely"]
    class_display = {
        "most likely": "Most Likely",
        "probable": "Probable",
        "least likely": "Least Likely",
    }
    class_palette = {
        "most likely": "#4C78A8",
        "probable":    "#F58518",
        "least likely": "#54A24B",
    }

    plot_df = selected_columns[top_features + ["label_encoded"]].copy()
    plot_df["label"] = plot_df["label_encoded"].map(label_map)

    n_feats = len(top_features)
    ncols = 3
    nrows = (n_feats + ncols - 1) // ncols

    def _add_sig_brackets(ax, group_data, n_comparisons=3):
        """Draw pairwise Mann-Whitney U significance brackets on ax.

        group_data : list of 1-D arrays, one per class in class_order order.
        n_comparisons : total comparisons per feature for Bonferroni correction.
        """
        all_vals = np.concatenate([d for d in group_data if len(d) > 0])
        all_vals = all_vals[np.isfinite(all_vals)]
        if len(all_vals) == 0:
            return
        y_min = float(np.nanmin(all_vals))
        y_max = float(np.nanmax(all_vals))
        y_range = float(y_max - y_min)
        if y_range <= 0:
            y_range = abs(y_max) if y_max != 0 else 1.0
        step = y_range * 0.12
        current_y = y_max + step * 0.35

        for i, j in [(0, 1), (0, 2), (1, 2)]:
            d1, d2 = group_data[i], group_data[j]
            if len(d1) < 3 or len(d2) < 3:
                current_y += step * 0.5
                continue
            _, p = _mannwhitneyu(d1, d2, alternative="two-sided")
            p_bonf = min(p * n_comparisons, 1.0)
            if p_bonf < 0.001:
                sig = "***"
            elif p_bonf < 0.01:
                sig = "**"
            elif p_bonf < 0.05:
                sig = "*"
            else:
                current_y += step * 0.5
                continue
            x1, x2 = float(i), float(j)
            bar_y = current_y + step * 0.1
            ax.plot([x1, x1, x2, x2], [current_y, bar_y, bar_y, current_y],
                    lw=0.8, color="black")
            ax.text((x1 + x2) / 2, bar_y, sig, ha="center", va="bottom", fontsize=8)
            current_y += step * 0.55
        top_limit = max(current_y + step * 0.35, y_max + step * 0.8)
        bottom_limit = y_min - step * 0.35
        ax.set_ylim(bottom=bottom_limit, top=top_limit)

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * 3.8))
    axes_flat = np.array(axes).flatten()

    for idx, feat in enumerate(top_features):
        ax = axes_flat[idx]
        group_data = []
        for xi, cls in enumerate(class_order):
            vals = plot_df.loc[plot_df["label"] == cls, feat].dropna().values
            group_data.append(vals)
            bp = ax.boxplot(
                vals,
                positions=[xi],
                widths=0.5,
                patch_artist=True,
                medianprops=dict(color="black", linewidth=1.5),
                whiskerprops=dict(linewidth=0.8),
                capprops=dict(linewidth=0.8),
                flierprops=dict(
                    marker="o", markersize=2, alpha=0.4,
                    markerfacecolor=class_palette[cls],
                    markeredgecolor="none",
                ),
            )
            for patch in bp["boxes"]:
                patch.set_facecolor(class_palette[cls])
                patch.set_alpha(0.75)

        ax.set_xticks(range(len(class_order)))
        ax.set_xticklabels([class_display[c] for c in class_order], rotation=20, ha="right", fontsize=8)
        ax.set_title(feat, fontsize=9, pad=4)
        ax.set_ylabel("Value", fontsize=7)
        ax.tick_params(axis="y", labelsize=7)
        _add_sig_brackets(ax, group_data, n_comparisons=3)
        if feat in {"Effect_SBP_Median", "PoPS_Score_SBP"}:
            y0, y1 = ax.get_ylim()
            ax.set_ylim(y0, y1 + (y1 - y0) * 0.12)

    for idx in range(n_feats, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    import matplotlib.patches as _mpatches
    legend_handles = [
        _mpatches.Patch(facecolor=class_palette[cls], alpha=0.75, label=class_display[cls])
        for cls in class_order
    ]
    fig.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(1.0, 0.0),
               ncol=1, fontsize=9, title="Gene class", title_fontsize=9)

    fig.suptitle(
        "Feature Distributions by Gene Class\n"
        "(features ranked by ANOVA F-statistic; "
        "* p<0.05  ** p<0.01  *** p<0.001, Bonferroni-corrected pairwise Mann-Whitney U)",
        y=1.01, fontsize=10,
    )
    plt.tight_layout()
    dist_plot_path = eda_out("top_feature_distributions_by_class.png")
    plt.savefig(dist_plot_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved feature distribution plot to: {dist_plot_path}")

# Publication plot 8: variance diagnostics on retained features
feature_var = X_retained.var(axis=0)
var_df = feature_var.reset_index()
var_df.columns = ["feature", "variance"]
var_df = var_df.sort_values("variance")
var_df.to_csv(eda_out("retained_feature_variance.tsv"), sep="\t", index=False)

nzv_threshold = 1e-4
nzv_df = var_df[var_df["variance"] <= nzv_threshold]
nzv_df.to_csv(eda_out("near_zero_variance_features.tsv"), sep="\t", index=False)

plt.figure(figsize=(8, 5))
plt.hist(var_df["variance"], bins=30)
plt.axvline(nzv_threshold, linestyle="--", color="red", linewidth=1)
plt.xlabel("Feature variance")
plt.ylabel("Number of features")
plt.title("Retained Feature Variance Distribution")
plt.tight_layout()
var_plot_path = eda_out("retained_feature_variance_hist.png")
plt.savefig(var_plot_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved variance histogram to: {var_plot_path}")
print(f"Saved near-zero variance feature list to: {eda_out('near_zero_variance_features.tsv')}")

# VIF calculation on retained features (multicollinearity diagnostic).
# A constant column is required by statsmodels' VIF implementation.
X_vif = X_retained.copy().dropna()
if X_vif.shape[0] > 0 and X_vif.shape[1] >= 2:
    X_vif_arr = X_vif.values.astype(float)
    # Prepend intercept column expected by variance_inflation_factor
    X_vif_with_const = np.column_stack([np.ones(X_vif_arr.shape[0]), X_vif_arr])
    vif_values = [
        variance_inflation_factor(X_vif_with_const, i + 1)
        for i in range(X_vif_arr.shape[1])
    ]
    vif_df = pd.DataFrame({"feature": X_vif.columns, "VIF": vif_values})
    vif_df = vif_df.sort_values("VIF", ascending=False)
    vif_path = eda_out("retained_feature_vif.tsv")
    vif_df.to_csv(vif_path, sep="\t", index=False)
    print(f"Saved VIF table to: {vif_path}")

    # Bar chart: colour by VIF > 10 (red), > 5 (orange), <= 5 (blue)
    top_vif_n = min(30, len(vif_df))
    top_vif = vif_df.head(top_vif_n).iloc[::-1]
    colors = [
        "#E45756" if v > 10 else "#F58518" if v > 5 else "#4C78A8"
        for v in top_vif["VIF"]
    ]
    plt.figure(figsize=(10, max(6, top_vif_n * 0.32)))
    plt.barh(top_vif["feature"], top_vif["VIF"], color=colors)
    plt.axvline(5, linestyle="--", color="#F58518", linewidth=1, label="VIF = 5")
    plt.axvline(10, linestyle="--", color="#E45756", linewidth=1, label="VIF = 10")
    plt.xlabel("Variance Inflation Factor (VIF)")
    plt.ylabel("Feature")
    plt.title("VIF of Retained Features (top 30)")
    plt.legend()
    plt.tight_layout()
    vif_plot_path = eda_out("retained_feature_vif.png")
    plt.savefig(vif_plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved VIF plot to: {vif_plot_path}")
else:
    print("Skipping VIF: insufficient data after dropping NaN rows.")

ensure_parent_dir(config.cleaned_training_data_eda)
ensure_parent_dir(config.cleaned_training_data_ml)
selected_columns.to_csv(config.cleaned_training_data_eda)
selected_columns.to_csv(config.cleaned_training_data_ml)
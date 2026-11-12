"""K-BRIDGE and the within-study comparator classifiers."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

# Preserve the original floating-point handling for the BLAS back end.
np.seterr(all="ignore")
DATA = Path(__file__).with_name("data")

SEED = 42


MODULES = {
    "BCR-proximal": [
        "CD79A", "CD79B", "BLNK", "BTK", "SYK", "LYN", "FYN", "PLCG2",
        "PIK3CD", "PIK3AP1", "CARD11", "MALT1", "BCL10",
    ],
    "Immediate-early": [
        "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3",
        "NR4A1", "DUSP1", "CD69", "MYC",
    ],
    "NF-kB": [
        "NFKB1", "RELA", "REL", "NFKBIA", "NFKBIZ", "TNFAIP3", "TRAF6",
        "IKBKB", "IKBKG",
    ],
    "Memory-B": [
        "CD27", "MS4A1", "CD37", "CD44", "CD82", "BANK1", "TNFRSF13B",
        "CD22", "CR2", "PAX5", "BACH2", "SPIB",
    ],
    "ASC-differentiation": [
        "PRDM1", "XBP1", "IRF4", "SDC1", "MZB1", "SEC11C", "DERL3",
        "JCHAIN", "TNFRSF17", "SLAMF7", "ELL2",
    ],
    "Antibody-secretion": [
        "IGKC", "IGLC1", "IGLC2", "IGHA1", "IGHG1", "IGHM", "IGLL5",
        "SSR4", "SEC61A1", "HSPA5",
    ],
    "B-cell-identity": [
        "CD19", "CD79A", "CD79B", "MS4A1", "CD22", "CD74", "HLA-DRA",
        "CD37", "BANK1", "PAX5", "SPIB",
    ],
}


def bootstrap_metric(y: np.ndarray, p: np.ndarray, metric, n_boot: int = 1000) -> tuple[float, float, float]:
    """Class-stratified nonparametric bootstrap for a binary metric."""
    vals = []
    rng = np.random.default_rng(SEED)
    class_indices = [np.flatnonzero(y == label) for label in np.unique(y)]
    for _ in range(n_boot):
        ix = np.concatenate([
            rng.choice(indices, size=len(indices), replace=True)
            for indices in class_indices
        ])
        vals.append(metric(y[ix], p[ix]))
    return float(metric(y, p)), float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


def module_incidence(genes: list[str]) -> tuple[np.ndarray, list[str]]:
    names = list(MODULES)
    mat = np.zeros((len(genes), len(names)))
    for j, name in enumerate(names):
        members = set(MODULES[name])
        for i, gene in enumerate(genes):
            mat[i, j] = gene in members
    return mat, names


def normalized_graph(genes: list[str], use_graph: bool = True, alpha: float = 0.35) -> np.ndarray:
    if not use_graph:
        return np.eye(len(genes))
    m, _ = module_incidence(genes)
    a = m @ m.T
    np.fill_diagonal(a, 0)
    a = a + np.eye(len(genes))
    degree = np.maximum(a.sum(axis=1), 1)
    norm = a / np.sqrt(np.outer(degree, degree))
    return (1 - alpha) * np.eye(len(genes)) + alpha * norm


class KBridge:
    def __init__(self, genes: list[str], use_graph=True, use_gate=True,
                 gene_features=True, pathway_features=True, seed=SEED):
        self.genes = genes
        self.use_graph = use_graph
        self.use_gate = use_gate
        self.gene_features = gene_features
        self.pathway_features = pathway_features
        self.seed = seed
        self.S = normalized_graph(genes, use_graph=use_graph)
        self.M, self.module_names = module_incidence(genes)
        self.Mn = self.M / np.maximum(self.M.sum(axis=0, keepdims=True), 1)

    def fit(self, x: np.ndarray, y: np.ndarray):
        self.scaler = StandardScaler().fit(x)
        z = self.scaler.transform(x)
        mi = mutual_info_classif(z, y, random_state=self.seed)
        if np.max(mi) > np.min(mi):
            mi = (mi - np.min(mi)) / (np.max(mi) - np.min(mi))
        else:
            mi = np.zeros_like(mi)
        self.gate = 0.5 + mi if self.use_gate else np.ones_like(mi)
        f = self._features_from_z(z)
        self.clf = LogisticRegression(
            penalty="elasticnet", solver="saga", l1_ratio=0.5, C=0.5,
            class_weight="balanced", max_iter=8000, random_state=self.seed,
        ).fit(f, y)
        return self

    def _features_from_z(self, z: np.ndarray) -> np.ndarray:
        # Some BLAS builds raise floating-point status flags for finite products.
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            d = (z @ self.S) * self.gate
            blocks = []
            if self.gene_features:
                blocks.append(d)
            if self.pathway_features:
                blocks.append(d @ self.Mn)
        features = np.column_stack(blocks)
        if not np.isfinite(features).all():
            raise ValueError("K-BRIDGE produced non-finite features")
        return features

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        return self.clf.predict_proba(self._features_from_z(self.scaler.transform(x)))

    def project_z(self, z: np.ndarray) -> np.ndarray:
        """Cross-domain relative projection from within-cohort standardized genes."""
        return self.clf.decision_function(self._features_from_z(z))

    def feature_names(self) -> list[str]:
        names = []
        if self.gene_features:
            names.extend(self.genes)
        if self.pathway_features:
            names.extend([f"module:{x}" for x in self.module_names])
        return names


def load_human_training() -> tuple[pd.DataFrame, pd.Series]:
    metadata = pd.read_csv(DATA / "GSE214265_metadata.csv.gz", sep=";", decimal=",", index_col=0)
    expression = pd.read_csv(DATA / "GSE214265_expression.txt.gz", sep="\t", index_col=0,
                             encoding="latin-1")
    expression.index = expression.index.astype(str).str.upper()
    expression = expression.groupby(level=0).mean()
    cells = [c for c in metadata.index if c in expression.columns]
    genes = sorted(set(sum(MODULES.values(), [])) & set(expression.index))
    x = np.log1p(expression.loc[genes, cells].T.astype(float))
    y = (metadata.loc[cells, "state"] == "ASC").astype(int)
    return x, y


def evaluate_models(x: pd.DataFrame, y: pd.Series):
    genes = list(x.columns)
    splits = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(x, y))
    models = {
        "K-BRIDGE": lambda seed: KBridge(genes, seed=seed),
        "Elastic-net": lambda seed: make_pipeline(
            StandardScaler(), LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5,
                                                  C=0.5, class_weight="balanced", max_iter=8000,
                                                  random_state=seed)),
        "RBF-SVM": lambda seed: make_pipeline(StandardScaler(), SVC(C=2, kernel="rbf", probability=True,
                                                                     class_weight="balanced", random_state=seed)),
        "Random forest": lambda seed: RandomForestClassifier(n_estimators=500, min_samples_leaf=3,
                                                               class_weight="balanced", random_state=seed,
                                                               n_jobs=-1),
        "XGBoost 2.1": lambda seed: XGBClassifier(n_estimators=300, max_depth=3, learning_rate=0.03,
                                                   subsample=0.8, colsample_bytree=0.8,
                                                   eval_metric="logloss", random_state=seed,
                                                   n_jobs=4),
        "MLP": lambda seed: make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(32, 16),
                                                                            alpha=0.01, max_iter=1500,
                                                                            early_stopping=True,
                                                                            random_state=seed)),
    }
    predictions = {}
    rows = []
    for name, make_model in models.items():
        oof = np.zeros(len(y))
        for fold, (tr, te) in enumerate(splits):
            model = make_model(SEED + fold)
            model.fit(x.iloc[tr].values, y.iloc[tr].values)
            oof[te] = model.predict_proba(x.iloc[te].values)[:, 1]
        predictions[name] = oof
        auroc, roc_lo, roc_hi = bootstrap_metric(y.values, oof, roc_auc_score)
        auprc, pr_lo, pr_hi = bootstrap_metric(y.values, oof, average_precision_score)
        pred = (oof >= 0.5).astype(int)
        rows.append({
            "model": name, "AUROC": auroc, "AUROC_low": roc_lo, "AUROC_high": roc_hi,
            "AUPRC": auprc, "AUPRC_low": pr_lo, "AUPRC_high": pr_hi,
            "balanced_accuracy": balanced_accuracy_score(y, pred), "F1": f1_score(y, pred),
        })
    metrics = pd.DataFrame(rows).sort_values("AUROC", ascending=False)
    oof_df = pd.DataFrame(predictions, index=x.index)
    oof_df["label"] = y
    return metrics, oof_df, splits


def evaluate_ablation(x, y, splits):
    genes = list(x.columns)
    variants = {
        "Full": dict(),
        "No graph": dict(use_graph=False),
        "No gate": dict(use_gate=False),
        "Gene only": dict(pathway_features=False),
        "Pathway only": dict(gene_features=False),
    }
    rows = []
    for name, kwargs in variants.items():
        oof = np.zeros(len(y))
        for fold, (tr, te) in enumerate(splits):
            m = KBridge(genes, seed=SEED + fold, **kwargs).fit(x.iloc[tr].values, y.iloc[tr].values)
            oof[te] = m.predict_proba(x.iloc[te].values)[:, 1]
        rows.append({"variant": name, "AUROC": roc_auc_score(y, oof),
                     "AUPRC": average_precision_score(y, oof)})
    return pd.DataFrame(rows)


def seed_and_permutation_tests(x, y):
    genes = list(x.columns)
    seed_rows = []
    for seed in range(10):
        oof = np.zeros(len(y))
        folds = StratifiedKFold(5, shuffle=True, random_state=100 + seed)
        for fold, (tr, te) in enumerate(folds.split(x, y)):
            m = KBridge(genes, seed=1000 + seed * 10 + fold).fit(x.iloc[tr].values, y.iloc[tr].values)
            oof[te] = m.predict_proba(x.iloc[te].values)[:, 1]
        seed_rows.append({"seed": seed, "AUROC": roc_auc_score(y, oof),
                          "AUPRC": average_precision_score(y, oof)})
    null = []
    fixed_splits = list(StratifiedKFold(5, shuffle=True, random_state=SEED).split(x, y))
    rng = np.random.default_rng(SEED)
    for perm in range(100):
        yp = pd.Series(rng.permutation(y.values), index=y.index)
        oof = np.zeros(len(y))
        for fold, (tr, te) in enumerate(fixed_splits):
            m = KBridge(genes, seed=2000 + perm * 10 + fold).fit(x.iloc[tr].values, yp.iloc[tr].values)
            oof[te] = m.predict_proba(x.iloc[te].values)[:, 1]
        null.append(roc_auc_score(yp, oof))
    return pd.DataFrame(seed_rows), np.asarray(null)

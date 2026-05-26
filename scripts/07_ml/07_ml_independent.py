import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.metrics import RocCurveDisplay, accuracy_score, f1_score, roc_auc_score
import xgboost as xgb
import shap

TABLE_DIR = Path("results/tables")
FIG_DIR   = Path("results/figures")
LOG_DIR   = Path("results/logs")

def main():
    print("="*55)
    print("STEP 07 - ML AMR Prediction (Independent Labels)")
    print(f"Date: {datetime.now().isoformat()}")
    print("="*55)

    # Load resistance matrix
    df = pd.read_csv(TABLE_DIR / "resistome_matrix.tsv", sep="\t", index_col=0)
    print(f"Matrix shape: {df.shape}")
    print(f"Columns: {df.columns.tolist()}")

    # Identify numeric feature columns (exclude non-numeric metadata)
    feature_cols = []
    for col in df.columns:
        try:
            # Check if column can be converted to numeric
            pd.to_numeric(df[col])
            feature_cols.append(col)
        except (ValueError, TypeError):
            print(f"Excluding non-numeric column: {col}")
    
    X = df[feature_cols].values.astype(float)
    print(f"Using {len(feature_cols)} feature columns")

    # Load independent labels
    labels = pd.read_csv(TABLE_DIR / "amrfinder_mdr_labels_independent.tsv", sep="\t")
    print(f"Labels shape: {labels.shape}")

    # Align indices – assume df index contains accessions
    # If df index is not accession, you may need to adjust
    common = df.index.intersection(labels["accession"])
    X = X[df.index.isin(common)]
    y = labels.set_index("accession").loc[common]["MDR_amrfinder"].values

    print(f"Samples:  {len(y)}")
    print(f"MDR (independent): {y.sum()} ({y.mean()*100:.1f}%)")
    print(f"Non-MDR:  {len(y)-y.sum()} ({(1-y.mean())*100:.1f}%)")

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    scoring = ["accuracy", "f1_weighted", "roc_auc"]

    print("\n--- Random Forest ---")
    rf = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
    rf_scores = cross_validate(rf, X, y, cv=skf, scoring=scoring, error_score='raise')
    print(f"  Accuracy: {rf_scores['test_accuracy'].mean():.3f} +/- {rf_scores['test_accuracy'].std():.3f}")
    print(f"  F1:       {rf_scores['test_f1_weighted'].mean():.3f} +/- {rf_scores['test_f1_weighted'].std():.3f}")
    print(f"  ROC-AUC:  {rf_scores['test_roc_auc'].mean():.3f} +/- {rf_scores['test_roc_auc'].std():.3f}")

    print("\n--- XGBoost ---")
    xgb_model = xgb.XGBClassifier(n_estimators=200, random_state=42,
                                    eval_metric="logloss", verbosity=0)
    xgb_scores = cross_validate(xgb_model, X, y, cv=skf, scoring=scoring, error_score='raise')
    print(f"  Accuracy: {xgb_scores['test_accuracy'].mean():.3f} +/- {xgb_scores['test_accuracy'].std():.3f}")
    print(f"  F1:       {xgb_scores['test_f1_weighted'].mean():.3f} +/- {xgb_scores['test_f1_weighted'].std():.3f}")
    print(f"  ROC-AUC:  {xgb_scores['test_roc_auc'].mean():.3f} +/- {xgb_scores['test_roc_auc'].std():.3f}")

    # Fit full models for SHAP and ROC curves
    print("\nFitting full models for SHAP and ROC curves...")
    rf.fit(X, y)
    xgb_model.fit(X, y)

    # Plot ROC curves
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Figure 3: ML model performance - MDR prediction (Independent Labels)", fontsize=13, fontweight="bold")

    RocCurveDisplay.from_estimator(rf, X, y, ax=axes[0], name="Random Forest", color="#3C3489")
    axes[0].plot([0,1],[0,1],"k--",alpha=0.4)
    axes[0].set_title("A. Random Forest ROC curve", fontsize=11, fontweight="bold")

    RocCurveDisplay.from_estimator(xgb_model, X, y, ax=axes[1], name="XGBoost", color="#1D9E75")
    axes[1].plot([0,1],[0,1],"k--",alpha=0.4)
    axes[1].set_title("B. XGBoost ROC curve", fontsize=11, fontweight="bold")

    plt.tight_layout()
    plt.savefig(FIG_DIR / "figure3_roc_curves_independent.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIG_DIR / "figure3_roc_curves_independent.pdf", bbox_inches="tight")
    plt.close()
    print("Figure 3 saved")

    # SHAP analysis
    print("Calculating SHAP values...")
    explainer = shap.TreeExplainer(rf)
    shap_values = explainer.shap_values(X)
    if isinstance(shap_values, list):
        shap_values = shap_values[1]  # for binary classification

    plt.figure(figsize=(10, 6))
    shap.summary_plot(shap_values, X, feature_names=feature_cols,
                      plot_type="bar", max_display=15, show=False)
    plt.title("Figure 4: SHAP feature importance - Random Forest (Independent Labels)", fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "figure4_shap_importance_independent.png", dpi=300, bbox_inches="tight")
    plt.savefig(FIG_DIR / "figure4_shap_importance_independent.pdf", bbox_inches="tight")
    plt.close()
    print("Figure 4 saved")

    # Save results
    summary = {
        "date": datetime.now().isoformat(),
        "n_genomes": len(y),
        "n_features": len(feature_cols),
        "mdr_count": int(y.sum()),
        "mdr_pct": round(float(y.mean()*100),1),
        "independent_labels_source": "NCBI AMRFinderPlus",
        "random_forest": {
            "accuracy": round(float(rf_scores["test_accuracy"].mean()),3),
            "f1": round(float(rf_scores["test_f1_weighted"].mean()),3),
            "roc_auc": round(float(rf_scores["test_roc_auc"].mean()),3),
        },
        "xgboost": {
            "accuracy": round(float(xgb_scores["test_accuracy"].mean()),3),
            "f1": round(float(xgb_scores["test_f1_weighted"].mean()),3),
            "roc_auc": round(float(xgb_scores["test_roc_auc"].mean()),3),
        },
    }
    with open(LOG_DIR / "07_ml_results_independent.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nLog saved: {LOG_DIR / '07_ml_results_independent.json'}")
    print("\n=== MANUSCRIPT TABLE 3 (Independent Labels) ===")
    print(f"{'Model':<20} {'Accuracy':>10} {'F1':>10} {'ROC-AUC':>10}")
    print("-"*52)
    print(f"{'Random Forest':<20} {rf_scores['test_accuracy'].mean():>10.3f} {rf_scores['test_f1_weighted'].mean():>10.3f} {rf_scores['test_roc_auc'].mean():>10.3f}")
    print(f"{'XGBoost':<20} {xgb_scores['test_accuracy'].mean():>10.3f} {xgb_scores['test_f1_weighted'].mean():>10.3f} {xgb_scores['test_roc_auc'].mean():>10.3f}")

if __name__ == "__main__":
    main()
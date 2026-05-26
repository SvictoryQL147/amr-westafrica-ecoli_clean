# AMR West Africa Genomics Study

## Genomic Characterization and Machine Learning‑Based Prediction of Antimicrobial Resistance in *Escherichia coli* from West Africa

**Author:** Chukwuka Onyekwere  
**Institution:** Nigeria  
**Start Date:** 2026-05-23  
**Last Updated:** 2026-05-26  
**Target Journal:** Microbial Genomics / BMC Genomics  
**Status:** ✅ Complete — manuscript ready for submission

---

## Study Overview

This study performs whole‑genome sequence (WGS) analysis of 235 *Escherichia coli* isolates from Nigeria (n=185) and Ghana (n=50) sourced from NCBI public databases. It characterises the resistome landscape, detects clinically critical resistance genes (including blaNDM‑1), and develops machine learning models to predict multidrug resistance (MDR) using independent AMRFinderPlus gold standard labels.

### Research Questions (Answered)

1. **What is the resistome landscape of West African *E. coli* isolates?**  
   → Dominated by blaTEM (55.3%), sul2 (48.5%), qnrS (26.8%). ESBL blaCTX‑M (4.3%), carbapenem blaNDM‑1 (0.9%) detected.

2. **What is the MDR prevalence and are there country‑level differences?**  
   → Independent AMRFinderPlus labels: 52.8% MDR overall; Ghana 64%, Nigeria 50%. K‑mer screening gave 22.6% (circular). Ghanaian isolates show higher rates, possibly due to cohort effects.

3. **Can ML models accurately predict MDR phenotype from genotype without circularity?**  
   → Yes. Random Forest and XGBoost predict independent AMRFinderPlus MDR with ROC‑AUC = 0.954 (±0.011). SHAP identifies blaTEM, sul2, qnrS as top predictors.

### Hypothesis (Revised)

West African *E. coli* isolates harbor a distinct resistome profile with high prevalence of beta‑lactam, sulfonamide, and fluoroquinolone resistance determinants. MDR can be reliably predicted from a compact gene panel using independent gold standard labels.

---

## Project Structure
amr_westAfrica_study/
├── data/
│ ├── raw/ # Downloaded genome files (FASTA) – not included in repo
│ ├── processed/ # TSV files from NCBI Pathogen Detection
│ └── metadata/ # master_metadata.tsv, gca_to_gcf.json
├── scripts/
│ ├── 02_metadata/ # Metadata parsing
│ ├── 04_resistome/ # Resistome profiling & independent label generation
│ └── 07_ml/ # Machine learning (independent labels)
├── results/
│ ├── tables/ # resistome_matrix.tsv, amrfinder_mdr_labels_independent.tsv
│ ├── figures/ # ROC curves, SHAP plots
│ └── logs/ # Analysis logs
├── manuscript/ # Draft sections
│ ├── MANUSCRIPT_DRAFT.md
│ ├── abstract.md
│ ├── introduction.md
│ └── discussion.md
├── LAB_NOTEBOOK.md # Complete scientific diary
├── README.md # This file
└── .gitignore

text

---

## Pipeline Overview (Completed Steps)

| Step | Tool / Method | Input | Output | Status |
|------|---------------|-------|--------|--------|
| 1. Data collection | NCBI Datasets CLI | Search terms (Nigeria, Ghana) | 235 *E. coli* genomes (FASTA) | ✅ Complete |
| 2. Quality control | CheckM | Raw assemblies | Completeness ≥95%, contamination ≤5% | ✅ Complete |
| 3. Metadata parsing | Python (pandas) | JSON from NCBI | master_metadata.tsv | ✅ Complete |
| 4. Resistome profiling | k‑mer screening (21 genes) | FASTA | resistome_matrix.tsv | ✅ Complete |
| 5. Independent MDR labels | NCBI AMRFinderPlus (TSV) | BioSample IDs | amrfinder_mdr_labels_independent.tsv | ✅ Complete |
| 6. ML prediction | Random Forest, XGBoost, SHAP | 21‑gene matrix + independent labels | ROC‑AUC 0.954, SHAP importance | ✅ Complete |

**Note:** Pangenome, phylogenetic, and plasmid analyses were not performed (out of scope). The study focuses on resistome and ML prediction.

---

## Key Results Summary

| Finding | Result |
|---------|--------|
| Most prevalent genes (k‑mer) | blaTEM (55.3%), sul2 (48.5%), qnrS (26.8%) |
| ESBL (blaCTX‑M) | 4.3% |
| Carbapenemase (blaNDM‑1) | 2 Nigerian isolates (0.9%) |
| Independent MDR (AMRFinderPlus) | 124/235 (52.8%) |
| Random Forest ROC‑AUC | 0.954 (±0.011) |
| XGBoost ROC‑AUC | 0.954 (±0.009) |
| Top SHAP predictors | blaTEM, sul2, qnrS |

**Surveillance gap:** No *E. coli* WGS assemblies from Senegal, Cameroon, or Côte d'Ivoire in NCBI (as of May 2026).

---

## How to Reproduce This Study

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/amr-westafrica-ecoli.git
cd amr-westafrica-ecoli
2. Set up the environment (optional – scripts use standard Python libraries)
bash
# Using conda (if you have environment file)
conda env create -f envs/amr_env.yml
conda activate amr_study

# Or install packages with pip
pip install pandas numpy matplotlib scikit-learn xgboost shap
3. Run the independent MDR label parser
bash
python scripts/04_resistome/12_parse_tsv_by_biosample.py
4. Run ML training and evaluation
bash
python scripts/07_ml/07_ml_independent.py
All outputs will be written to results/tables/ and results/figures/.

Manuscript
The final manuscript is available in MANUSCRIPT_DRAFT.md. All sections (abstract, introduction, methods, results, discussion) are complete and corrected.

Changelog
Date	Entry
2026-05-23	Project initialized. Directory structure created.
2026-05-23	235 E. coli genomes downloaded (Nigeria 185, Ghana 50).
2026-05-23	Resistome profiling (k‑mer screening) completed.
2026-05-26	Circular ML (perfect AUC) identified and corrected.
2026-05-26	Independent AMRFinderPlus labels generated.
2026-05-26	Final ML with independent labels: ROC‑AUC 0.954.
2026-05-26	Manuscript updated, phylogenetics removed, ready for submission.
License
[Choose a license, e.g., MIT or CC‑BY‑4.0]

Contact
Chukwuka Onyekwere – [your.email@institution.edu.ng]
Collaborator – [friend's email]

For questions or collaboration, please open an issue on GitHub.

text

This README is now honest, up‑to‑date, and ready for your repository. It correctly reflects that the study is complete and uses *E. coli*, not *Klebsiella*. Push it along with the other corrected files.
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

TABLE_DIR = Path("results/tables")
FIG_DIR   = Path("results/figures")

df = pd.read_csv(TABLE_DIR / "amrfinder_raw.tsv", sep="\t")

gene_cols = ["blaTEM","sul2","qnrS","sul1","blaCTX-M","aph(3')-Ia",
             "blaNDM-1","tetA","tetB","dfrA1","aac(3)-II","qnrB",
             "blaOXA-1","mcr-1","blaSHV","rmtB","armA","blaCMY-2",
             "blaKPC-2","aac(6')-Ib"]

class_colors = {
    "blaTEM":"#3C3489","blaCTX-M":"#3C3489","blaSHV":"#3C3489",
    "blaOXA-1":"#3C3489","blaCMY-2":"#3C3489","blaKPC-2":"#3C3489","blaNDM-1":"#3C3489",
    "tetA":"#1D9E75","tetB":"#1D9E75",
    "sul1":"#EF9F27","sul2":"#EF9F27",
    "dfrA1":"#D85A30",
    "aac(3)-II":"#D4537E","aph(3')-Ia":"#D4537E","rmtB":"#D4537E",
    "armA":"#D4537E","aac(6')-Ib":"#D4537E",
    "qnrB":"#639922","qnrS":"#639922",
    "mcr-1":"#E24B4A",
}

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle("Figure 2: Resistome profile of E. coli from West Africa (n=235)",
             fontsize=13, fontweight="bold")

prev = df[gene_cols].mean().sort_values(ascending=True) * 100
prev = prev[prev > 0]
colors = [class_colors.get(g, "#888780") for g in prev.index]
bars = axes[0].barh(prev.index, prev.values, color=colors, edgecolor="white", height=0.7)
axes[0].set_xlabel("Prevalence (%)", fontsize=11)
axes[0].set_title("A. Resistance gene prevalence", fontsize=11, fontweight="bold")
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)
for bar, val in zip(bars, prev.values):
    axes[0].text(val + 0.5, bar.get_y() + bar.get_height()/2,
                 f"{val:.1f}%", va="center", fontsize=9)
axes[0].set_xlim(0, 70)

legend_elements = [
    mpatches.Patch(color="#3C3489", label="Beta-lactam"),
    mpatches.Patch(color="#1D9E75", label="Tetracycline"),
    mpatches.Patch(color="#EF9F27", label="Sulfonamide"),
    mpatches.Patch(color="#D85A30", label="Trimethoprim"),
    mpatches.Patch(color="#D4537E", label="Aminoglycoside"),
    mpatches.Patch(color="#639922", label="Fluoroquinolone"),
    mpatches.Patch(color="#E24B4A", label="Colistin"),
]
axes[0].legend(handles=legend_elements, loc="lower right", fontsize=8, title="Drug class")

countries = ["Nigeria", "Ghana", "Total"]
mdr_rates = [18.9, 36.0, 22.6]
non_rates = [100-r for r in mdr_rates]
x = np.arange(len(countries))
b1 = axes[1].bar(x, mdr_rates, color="#E24B4A", label="MDR", width=0.5, edgecolor="white")
b2 = axes[1].bar(x, non_rates, color="#B5D4F4", label="Non-MDR", width=0.5,
                 bottom=mdr_rates, edgecolor="white")
axes[1].set_xticks(x)
axes[1].set_xticklabels(countries, fontsize=11)
axes[1].set_ylabel("Percentage of genomes (%)", fontsize=11)
axes[1].set_title("B. MDR prevalence by country", fontsize=11, fontweight="bold")
axes[1].set_ylim(0, 110)
axes[1].legend(fontsize=10)
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
for bar, rate in zip(b1, mdr_rates):
    axes[1].text(bar.get_x() + bar.get_width()/2, rate/2,
                 f"{rate}%", ha="center", va="center",
                 color="white", fontsize=11, fontweight="bold")

plt.tight_layout()
plt.savefig(FIG_DIR / "figure2_resistome_profile.png", dpi=300, bbox_inches="tight")
plt.savefig(FIG_DIR / "figure2_resistome_profile.pdf", bbox_inches="tight")
plt.close()
print("Figure 2 saved:")
print("  results/figures/figure2_resistome_profile.png")
print("  results/figures/figure2_resistome_profile.pdf")

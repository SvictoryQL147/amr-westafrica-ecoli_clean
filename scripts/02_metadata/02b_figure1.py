import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("Figure 1: E. coli genome assemblies from West Africa",
             fontsize=13, fontweight="bold", y=1.02)

# Panel A: Genome counts by country
countries = ["Nigeria", "Ghana"]
counts    = [185, 50]
colors    = ["#3C3489", "#1D9E75"]
bars = axes[0].bar(countries, counts, color=colors, width=0.5, edgecolor="white")
axes[0].set_ylabel("Number of genome assemblies", fontsize=11)
axes[0].set_title("A. Genomes per country", fontsize=11, fontweight="bold")
axes[0].set_ylim(0, 220)
for bar, count in zip(bars, counts):
    axes[0].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 4,
                 str(count), ha="center", va="bottom", fontsize=11, fontweight="bold")
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)

# Panel B: Assembly level stacked bar
levels     = ["Scaffold", "Contig", "Chromosome"]
nigeria_v  = [118, 44, 23]
ghana_v    = [0, 50, 0]
level_cols = ["#7F77DD", "#1D9E75", "#EF9F27"]
x = np.arange(2)
bottom_n = np.zeros(2)
for level, col, nv, gv in zip(levels, level_cols, nigeria_v, ghana_v):
    vals = [nv, gv]
    axes[1].bar(x, vals, bottom=bottom_n, color=col,
                label=level, width=0.5, edgecolor="white")
    for j, v in enumerate(vals):
        if v > 0:
            axes[1].text(x[j], bottom_n[j] + v/2, str(v),
                         ha="center", va="center", fontsize=10,
                         color="white", fontweight="bold")
    bottom_n += np.array(vals)

axes[1].set_xticks(x)
axes[1].set_xticklabels(countries)
axes[1].set_ylabel("Number of genome assemblies", fontsize=11)
axes[1].set_title("B. Assembly level breakdown", fontsize=11, fontweight="bold")
axes[1].legend(title="Assembly level", bbox_to_anchor=(1.05, 1), loc="upper left")
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig("results/figures/figure1_genome_distribution.png", dpi=300, bbox_inches="tight")
plt.savefig("results/figures/figure1_genome_distribution.pdf", bbox_inches="tight")
plt.close()
print("Done.")
print("Saved: results/figures/figure1_genome_distribution.png (300 DPI)")
print("Saved: results/figures/figure1_genome_distribution.pdf (vector)")

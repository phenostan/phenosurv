"""Worked example: fit, diagnose and predict with ThreshSurvGDD.

Uses the simulated data in tests/data/ (4 site-years, 3 genotypes, 24
flowering events) and writes the README figures to docs/figures/.

Run it through pixi:

    pixi run example
"""

import json
from pathlib import Path

import matplotlib
import matplotlib.patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from phenosurv import ThreshSurvGDD

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data" / "testdata.json"
TRUTH = ROOT / "tests" / "data" / "testdata_truth.json"
FIG_DIR = ROOT / "docs" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# One fixed colour per genotype (never re-assigned), plus text/grid inks.
GENO_COLORS = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK_2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": INK_2,
    "axes.labelcolor": INK,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "font.size": 10,
    "legend.frameon": False,
})

# --- 1. Load data -------------------------------------------------------------
data = json.loads(DATA.read_text())
truth = json.loads(TRUTH.read_text())
events = np.array(data["events"])
precursors = np.array(data["precursor_events"])
genotype = np.array(data["genotype"])
site_year = np.repeat(np.arange(data["N_site_years"]), data["N_events"])

# --- 2. Fit and diagnose ------------------------------------------------------
model = ThreshSurvGDD()           # compiles the Stan model on first use
model.fit(data, seed=1)
print(model.diagnose())
summary = model.summary()
print(summary.loc[["T_base", "Psi0_bar", "tau_psi0", "sigma",
                   "Psi0[1]", "Psi0[2]", "Psi0[3]"]].round(2))

# --- 3. Predict ---------------------------------------------------------------
# (a) Retrodict the observed events.
pred_obs = model.predict(seed=2)                       # (draws, N)

# (b) New data: the same site-years and plantings under +2 C warming.
#     `events` is omitted, since it is what we are predicting.
warm = {k: v for k, v in data.items() if k != "events"}
warm["obs_temps"] = [t + 2.0 for t in data["obs_temps"]]
pred_warm = model.predict(warm, seed=3)

# --- Figure 1: the data ---------------------------------------------------------
# Per site-year: daily temperature on top, with the forcing it contributes
# between planting and the last flowering shaded (degrees above the true T_base). Below it, one row
# per plant: a bar from planting (precursor) to flowering (event).
fig = plt.figure(figsize=(10, 7))
outer = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.12)
temp_axes, strip_axes = [], []
for i in range(data["N_site_years"]):
    cell = outer[i // 2, i % 2].subgridspec(2, 1, height_ratios=[2.2, 1], hspace=0.08)
    ax = fig.add_subplot(cell[0], sharex=temp_axes[0] if temp_axes else None,
                         sharey=temp_axes[0] if temp_axes else None)
    strip = fig.add_subplot(cell[1], sharex=ax)
    temp_axes.append(ax)
    strip_axes.append(strip)

    start, n_days = data["temp_start_idxs"][i] - 1, data["N_days"][i]
    temps = np.array(data["obs_temps"][start:start + n_days])
    doy = np.arange(1, n_days + 1)
    rows = np.nonzero(site_year == i)[0]
    planting = precursors[rows].min()

    window = (doy >= planting) & (doy <= events[rows].max())
    ax.fill_between(doy, truth["T_base"], temps, where=window & (temps > truth["T_base"]),
                    color="#cde2fb", lw=0, interpolate=True)
    ax.plot(doy, temps, color="#8a8984", lw=1)
    ax.axhline(truth["T_base"], color=INK_2, lw=1, ls="--")
    ax.axvline(planting, color=INK, lw=1.5)
    ax.set_title(f"Site-year {i + 1}", loc="left")
    ax.tick_params(labelbottom=False)

    # Event strip: genotype rows, replicate plants offset within a row.
    rep = np.zeros(len(rows))
    for g in range(1, data["N_geno"] + 1):
        idx = np.nonzero(genotype[rows] == g)[0]
        rep[idx] = np.linspace(-0.18, 0.18, len(idx)) if len(idx) > 1 else 0
    for n, r in zip(rows, rep):
        y = genotype[n] + r
        color = GENO_COLORS[genotype[n] - 1]
        strip.plot([precursors[n], events[n]], [y, y], color=color, lw=1.5, alpha=0.5)
        strip.scatter(events[n], y, s=36, color=color, edgecolor="white",
                      linewidth=1.5, zorder=3)
    strip.axvline(planting, color=INK, lw=1.5)
    strip.set_ylim(data["N_geno"] + 0.6, 0.4)
    strip.set_yticks(range(1, data["N_geno"] + 1))
    strip.set_yticklabels([f"G{g}" for g in range(1, data["N_geno"] + 1)], fontsize=8)
    strip.grid(axis="y", visible=False)
    if i // 2 == 1:
        strip.set_xlabel("Day of year")
    if i % 2 == 0:
        ax.set_ylabel("Daily mean\ntemp. (°C)")
        strip.set_ylabel("Genotype")
    else:
        ax.tick_params(labelleft=False)
        strip.tick_params(labelleft=False)

temp_axes[0].set_xlim(1, 365)
handles = [plt.Line2D([], [], color=c, lw=1.5, marker="o", markersize=6,
                      markeredgecolor="white", label=f"Genotype {g + 1}: planting → flowering")
           for g, c in enumerate(GENO_COLORS)]
handles += [plt.Line2D([], [], color=INK, lw=1.5, label="Planting day"),
            plt.Line2D([], [], color=INK_2, lw=1, ls="--", label="True $T_{base}$ (5 °C)"),
            matplotlib.patches.Patch(color="#cde2fb", label="Forcing (degree-days) until last flowering")]
fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.995), fontsize=9)
fig.subplots_adjust(top=0.88, bottom=0.08, left=0.09, right=0.98)
fig.savefig(FIG_DIR / "data.png", dpi=150)

# --- Figure 2: posterior vs. the values the data were simulated from ----------
params = [("T_base", "$T_{base}$ (°C)", truth["T_base"], INK_2),
          ("sigma", "$\\sigma$ (degree-days)", truth["sigma"], INK_2)]
params += [(f"Psi0[{g + 1}]", f"$\\Psi_0$, genotype {g + 1} (degree-days)",
            truth["Psi0"][g], GENO_COLORS[g]) for g in range(data["N_geno"])]
fig, axes = plt.subplots(1, len(params), figsize=(13, 2.8))
draws = model.fit_obj.draws_pd()
for ax, (name, label, true, color) in zip(axes, params):
    ax.hist(draws[name], bins=40, color=color, alpha=0.85, edgecolor="white", linewidth=0.5)
    ax.axvline(true, color=INK, lw=2)
    ax.set_xlabel(label)
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["left"].set_visible(False)
axes[0].annotate("true value", (truth["T_base"], 1), xytext=(4, -10),
                 textcoords="offset points", xycoords=("data", "axes fraction"),
                 color=INK, fontsize=9)
fig.suptitle("Posterior distributions (histograms) recover the simulated values (black lines)",
             fontsize=11, fontweight="bold", x=0.01, ha="left")
fig.tight_layout()
fig.savefig(FIG_DIR / "posterior.png", dpi=150)

# --- Figure 3: retrodictive check ----------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 5))
lo, med, hi = np.percentile(pred_obs, [5, 50, 95], axis=0)
for g in range(data["N_geno"]):
    sel = genotype == g + 1
    ax.vlines(events[sel], lo[sel], hi[sel], color=GENO_COLORS[g], lw=2, alpha=0.6)
    ax.scatter(events[sel], med[sel], s=40, color=GENO_COLORS[g],
               edgecolor="white", linewidth=1.5, zorder=3, label=f"Genotype {g + 1}")
lims = [min(events.min(), lo.min()) - 5, max(events.max(), hi.max()) + 5]
ax.plot(lims, lims, color=INK_2, lw=1, ls="--")
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect("equal")
ax.set_xlabel("Observed flowering day")
ax.set_ylabel("Predicted flowering day\n(median and 90% interval)")
ax.set_title("Posterior retrodictive check", loc="left")
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(FIG_DIR / "retrodictive.png", dpi=150)

# --- Figure 4: prediction for new data (+2 C) --------------------------------
# One row per event: posterior median flowering day under the observed
# temperatures (grey) and the +2 C scenario (genotype colour).
med_now = np.median(pred_obs, axis=0)
med_warm = np.median(pred_warm, axis=0)
order = np.lexsort((genotype, site_year))[::-1]
fig, ax = plt.subplots(figsize=(7, 7.2))
for row, n in enumerate(order):
    color = GENO_COLORS[genotype[n] - 1]
    ax.plot([med_warm[n], med_now[n]], [row, row], color=GRID, lw=2, zorder=1)
    ax.scatter(med_now[n], row, s=40, color="#8a8984", edgecolor="white", linewidth=1.5, zorder=2)
    ax.scatter(med_warm[n], row, s=40, color=color, edgecolor="white", linewidth=1.5, zorder=3)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([f"Site-year {site_year[n] + 1} · G{genotype[n]}" for n in order], fontsize=8)
ax.grid(axis="y", visible=False)
ax.set_xlabel("Predicted flowering day (posterior median)")
ax.set_title("Predicted flowering: observed temps. (grey) vs. +2 °C (colour)", loc="left")
for g in range(data["N_geno"]):
    ax.scatter([], [], s=40, color=GENO_COLORS[g], label=f"Genotype {g + 1}, +2 °C")
ax.scatter([], [], s=40, color="#8a8984", label="Observed temperatures")
ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.09), ncol=2, fontsize=9)
fig.tight_layout()
fig.savefig(FIG_DIR / "prediction.png", dpi=150)

for g in range(data["N_geno"]):
    sel = genotype == g + 1
    print(f"Genotype {g + 1}: median shift under +2 C = "
          f"{np.median(med_warm[sel] - med_now[sel]):+.1f} days")
print(f"Figures written to {FIG_DIR}")

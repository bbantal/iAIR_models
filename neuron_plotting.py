#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 23:06:50 2023

@author: botond
"""


import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# raise
# rc params
plt.rcParams["font.size"] = 16
plt.rcParams["font.family"] = "Arial"

# Filepath
HOMEDIR = os.path.abspath(os.path.join(__file__, "../../../")) + "/"
SRCDIR = HOMEDIR + "data/axon_model/traub/"
OUTDIR = HOMEDIR + "results/axon_model/"

# Load data
# Spike trains
st1 = pd.read_csv(SRCDIR + "spiketrain_10.csv")
st2 = pd.read_csv(SRCDIR + "spiketrain_100.csv")

# Model results
df = pd.read_csv(SRCDIR + "traub_model_results.csv")

# Normalize to percentages
normalize = lambda series: (series-series.iloc[-1])/series.iloc[-1]*100*(-1 if series.iloc[-1]<0 else 1)
df = df.assign(**{"v_rest": normalize(df["v_rest"]),
             "freq": normalize(df["freq"]),
             "slope": normalize(df["slope"]),
             "v_peak": normalize(df["v_peak"])
             })

# Colors
colors = [None, None, "gold", "teal", "crimson", "deepskyblue"]

# Formatting props
s=1.5

# Gridspec
gridspec = mpl.gridspec.GridSpec(3, 2, height_ratios=(3, 1.5, 1.5))

# raise

# Create figure
plt.figure(figsize=(12, 10))

# %%
# Panel A
# plt.figure(figsize=(5, 4))
plt.subplot(gridspec[0])
plt.plot(st2["t"], st2["v"], color="gray", label=r"$v^{max}_{ATPase}=$" + "\n" + \
         "$100$ $µA/cm^2$", lw=0.5)
plt.plot(st1["t"], st1["v"], color="indianred", label=r"$v^{max}_{ATPase}=$" + "\n" + \
         "$10$ $µA/cm^2$", lw=0.5)
plt.xlim([0, 3500])
plt.ylim([-90, 35])
plt.xlabel("Time [ms]")
plt.ylabel("Membrane potential [mV]", labelpad=12)
plt.legend(fontsize=11, loc=1)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# Save
# plt.tight_layout()
# plt.savefig(OUTDIR + "traub_results_panelA.pdf") #, transparent=True)

# %%
# Panel B
# plt.figure(figsize=(5, 4))
plt.subplot(gridspec[1])
plt.plot(st1["t"], st1["v"], color="indianred", label="Ipumpmax=10")
plt.plot(st2["t"], st2["v"], color="gray", label="Ipumpmax=100")
plt.xlim([1500, 1550])
plt.ylim([-90, 35])
plt.xlabel("Time [ms]")
plt.ylabel("Membrane potential [mV]", labelpad=12)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# Save
# plt.tight_layout()
# plt.savefig(OUTDIR + "traub_results_panelB.pdf") #, transparent=True)

# %%
# Panel C
# plt.figure(figsize=(5, 3.2))
plt.subplot(gridspec[2])
sns.scatterplot(data=df, x="I_pump_max", y="v_rest",
                color=colors[2], linewidth=.5*s, edgecolor="black", s=100*s,
                zorder=2)
plt.plot(df["I_pump_max"], df["v_rest"], color=colors[2], lw=2*s, zorder=1)
plt.xlim([110, 1])
# plt.ylim([-88, -82])
plt.ylim([-1, 6])
plt.xlabel(r"$v^{max}_{ATPase}$" + " [" "$µA/cm^2$" + "]")
plt.ylabel("Relative change\nin resting membrane\npotential [%]", labelpad=16)
plt.xticks([100, 75, 50, 25])
plt.grid(zorder=0)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# Save
# plt.tight_layout()
# plt.savefig(OUTDIR + "traub_results_panelC.pdf") #, transparent=True)

# %%
# Panel D
# plt.figure(figsize=(5, 3.2))
plt.subplot(gridspec[3])
sns.scatterplot(data=df, x="I_pump_max", y="freq",
                color=colors[3], linewidth=.75*s, edgecolor="black", s=100*s, zorder=2)
plt.plot(df["I_pump_max"], df["freq"], color=colors[3], lw=2*s, zorder=1)
plt.xlim([110, 1])
# plt.ylim([55, 70])
plt.ylim([-5, 25])
plt.xlabel(r"$v^{max}_{ATPase}$" + " [" "$µA/cm^2$" + "]")
plt.ylabel("Relative change\nin firing frequency [%]", labelpad=20)
plt.xticks([100, 75, 50, 25])
plt.grid(zorder=0)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# Save
# plt.tight_layout()
# plt.savefig(OUTDIR + "traub_results_panelD.pdf") #, transparent=True)
# %%

# Panel E
# plt.figure(figsize=(5, 3.2))
plt.subplot(gridspec[4])
sns.scatterplot(data=df, x="I_pump_max", y="v_peak",
                color=colors[4], linewidth=.75*s, edgecolor="black", s=100*s, zorder=2)
plt.plot(df["I_pump_max"], df["v_peak"], color=colors[4], lw=2*s, zorder=1)
plt.xlim([110, 1])
# plt.ylim([27.4, 29.6])
# plt.yticks([27.5, 28.5, 29.5])
plt.ylim([-6.5, 1])
plt.yticks([-5, -2.5, 0])
plt.xlabel(r"$v^{max}_{ATPase}$" + " [" "$µA/cm^2$" + "]")
plt.ylabel("Relatove change\nin peak membrane\npotential [%]", labelpad=8)
plt.xticks([100, 75, 50, 25])
plt.grid(zorder=0)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# Save
# plt.tight_layout()
# plt.savefig(OUTDIR + "traub_results_panelE.pdf") #, transparent=True)
# %%
# Panel F
# plt.figure(figsize=(5, 3.2))
plt.subplot(gridspec[5])
sns.scatterplot(data=df, x="I_pump_max", y="slope",
                color=colors[5], linewidth=.75*s, edgecolor="black", s=100*s, zorder=2)
plt.plot(df["I_pump_max"], df["slope"], color=colors[5], lw=2*s, zorder=1)
plt.xlim([110, 1])
# plt.ylim([-0.21, -0.09])
plt.ylim([-110, 10])
plt.xlabel(r"$v^{max}_{ATPase}$" + " [" "$µA/cm^2$" + "]")
plt.ylabel("Relative change\nin slope of peak [%]", labelpad=2)
plt.xticks([100, 75, 50, 25])
plt.grid(zorder=0)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# Save
# plt.tight_layout()
# plt.savefig(OUTDIR + "traub_results_panelF.pdf") #, transparent=True)

# Tight layout
plt.tight_layout()

# Savefig
plt.savefig(OUTDIR + "traub_results.pdf", dpi=300)

# %%


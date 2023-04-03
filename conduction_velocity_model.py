#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 12:15:26 2023

@author: botond
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# rc params
plt.rcParams["font.size"] = 13
plt.rcParams["font.family"] = "Arial"

# Filepath
HOMEDIR = os.path.abspath(os.path.join(__file__, "../../../")) + "/"
SRCDIR = HOMEDIR + "data/axon_model/"
OUTDIR = HOMEDIR + "results/axon_model/"

# %%
# Compute data

# Model function
def calc_cv(V_r, V_a, V_t):
    cv = 1/(np.log((V_r - V_a)/(2*V_t - V_r - V_a)))
    return cv

# %%
# 1D figure - panels
# ------

# Generate data
# ----

V_t = -50.;
vr_def = -65.;
va_def = 30.;

# Compute reference value for nonscaled CV
ref_CV = calc_cv(vr_def, va_def, V_t)

# Grid
vr = -57 - 30*(np.logspace(np.log(1), np.log(11), 10, base=np.e)/10);
va = np.linspace(00, 50, 11);

# Compute non-scaled CV
CV_vr = calc_cv(vr, va_def, V_t)
CV_va = calc_cv(vr_def, va, V_t)

# Correct with reference
CV_corr_vr = (CV_vr - ref_CV)/ref_CV * 100;
CV_corr_va = (CV_va - ref_CV)/ref_CV * 100;

# Plot
# ------

colors = ["crimson", "deepskyblue"]
s=1.5

# plt.figure(figsize=(12, 4))
# plt.figure(figsize=(7.25, 4))

# Panel 1
# ---------

plt.figure(figsize=(5, 4))
# plt.subplot(1, 2, 1)
plt.scatter(va, CV_corr_va, color=colors[0], linewidth=.75*s,
            edgecolor="black", s=100*s, zorder=3)
plt.plot(va, CV_corr_va, color=colors[0], lw=2*s, zorder=2)

plt.xticks(np.arange(0, 55, 10))
plt.ylim([-75, 60])
plt.xlabel("Peak membrane\npotential [mV]")#, fontsize=13)
plt.ylabel("Relative change\nin conduction velocity (%)")#, fontsize=13)
plt.grid(zorder=1)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# # Save
plt.tight_layout()
plt.savefig(OUTDIR + "cv_formula_panel1.pdf") # transparent=True)

# Panel 2
# --------

plt.figure(figsize=(5, 4))
# plt.subplot(1, 2, 2)
plt.scatter(vr, CV_corr_vr, color=colors[1], linewidth=.75*s,
            edgecolor="black", s=100*s, zorder=3)
plt.plot(vr, CV_corr_vr, color=colors[1], lw=2*s, zorder=2)


plt.xticks(np.arange(-90, -55, 10))
plt.grid(zorder=1)
plt.xlim([-94, -56])
plt.ylim([-75, 60])
plt.xlabel("Resting membrane\npotential [mV]")#, fontsize=13)
plt.ylabel("Relative change\nin conduction velocity (%)")#, fontsize=13)

ax = plt.gca()
for sp in ['bottom', 'top', 'right', 'left']:
    ax.spines[sp].set_linewidth(1.4)
    ax.spines[sp].set_color("black")

# # Save
plt.tight_layout()
plt.savefig(OUTDIR + "cv_formula_panel2.pdf") #, transparent=True)


# plt.tight_layout(rect=[0.03, 0, .97, 1])
# plt.gcf().subplots_adjust(wspace=0.4)

# Save
# plt.savefig(OUTDIR + "cv_formula_1d.pdf", transparent=True)

#############





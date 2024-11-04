#!/usr/bin/env python
"""
Plots to compare PhaseTracer results against analytic expressions
"""

import os
import pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rc, cm


try:
    latex = os.environ["MATPLOTLIB_LATEX"]
except KeyError:
    latex = False
if latex:
    rc('text', usetex=True)
    rc('font', **{'family': 'serif', 'size': 14})

rc('axes', **{'grid': True})
rc('grid', **{'ls': ':'})


# Fetch data

data_file_name = 'Z2ScalarSingletModel_Results.txt'
data = pandas.read_csv(data_file_name, delim_whitespace=True)
failed = data.isnull().any(1)

l_hs = data["l_hs"]
m_s = data["m_s"]
TC = data["T1_c^PT"]
TC_theory = data["T_c^EX"]

print("{} failures in {} points".format(sum(failed), len(data["m_s"])))

# Plot lambda_hs against m_s

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6.4*2, 4.8))
cp = ax[0].scatter(m_s, l_hs, c=TC, s=2, edgecolor='None', cmap=cm.get_cmap('rainbow', 15), vmin=0, vmax=140, alpha=0.8)
# Add colorbar
cb = fig.colorbar(cp, ax=ax[0])
cb.set_label(r"$T_C$ (GeV)")


diff = abs(TC - TC_theory) / TC_theory
Z = zip(m_s, l_hs, diff)
Z.sort(key=lambda t: t[2], reverse=False)
m_s, l_hs, diff = zip(*Z)

cp = ax[1].scatter(m_s, l_hs, c=diff, marker="s", s=2, edgecolor='None', cmap=cm.get_cmap('autumn_r', 15))
# Add colorbar
cb = fig.colorbar(cp, ax=ax[1])
cb.set_label(r"Relative difference versus analytic")

for a in ax:
    a.set_ylabel(r"$\lambda_{hs}$")
    a.set_xlabel("$m_s$ (GeV)")
    a.set_xlim(10, 235)
    a.set_ylim(0.2, 2)


plt.savefig('ms_lambda_Z2_SSM.pdf', bbox_inches='tight')

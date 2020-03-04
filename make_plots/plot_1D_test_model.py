#!/usr/bin/env python
"""
Plots for 1D test model
=======================
"""


import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm, rc


try:
    latex = os.environ["MATPLOTLIB_LATEX"]
except KeyError:
    latex = False
if latex:
    rc('text', usetex=True)
    rc('font', **{'family': 'serif'})
    rc('font', **{'family': 'serif', 'size': 14})


def V(x, T):
    return -100*x**2 - 10.*x**3 + 0.1*x**4 + 0.1*x**2 * T**2


scale = 1e-6
TC = 59.2297
x = np.linspace(-20, 100, 600)
T = np.array([0, 20, 30, 40, 50, TC, 65, 70, 80])
xx, TT = np.meshgrid(x, T)
VV = V(xx, TT)

fig, ax = plt.subplots()
cp = ax.scatter(xx, scale * VV, c=TT, s=5, edgecolor='None', cmap=cm.get_cmap('rainbow', 15), alpha=0.8)

ax.text(43, -0.1, r'$T_C=59.2$ GeV')

plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 3))
plt.grid(True, linestyle=':')

ax.set_ylim(-1, 0.5)
ax.set_xlim(-20, 90)
ax.set_ylabel(r'$V(\phi) \times 10^{-6}$ (GeV)${}^4$')
ax.set_xlabel(r'$\phi$ (GeV)')

cb = fig.colorbar(cp)
cb.set_label(r"$T$ (GeV)")


plt.savefig('potential_1D_test_model.pdf', bbox_inches='tight')

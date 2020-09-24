#!/usr/bin/env python
"""
Compare h-bar expansion against https://arxiv.org/pdf/1808.01098.pdf fig. 2 
"""

from subprocess import run
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import subprocess


def call_pt(lambda_hs, Q, xi, tree_level_tadpoles):
    r = subprocess.run(["../bin/h_bar_expansion",
                        str(lambda_hs),
                        str(Q),
                        str(xi),
                        str(int(tree_level_tadpoles))],
                        stdout=subprocess.PIPE, encoding='utf8')
    TC = np.nan
    gamma = np.nan
    mh_tree = np.nan
    mh_1l = np.nan

    for line in r.stdout.split("\n"):
        if line.startswith("TC = "):
            TC = float(line.split("=")[-1].strip())
        if line.startswith("gamma_HT = "):
            gamma = float(line.split("=")[-1].strip())
        if line.startswith("mh_tree = "):
            mh_tree = float(line.split("=")[-1].strip())
        if line.startswith("mh_1l = "):
            mh_1l = float(line.split("=")[-1].strip())

    return TC, gamma, mh_tree, mh_1l

if __name__ == "__main__":

    fig, ax = plt.subplots(1, 3, figsize=(12, 4))

    lambda_hs = np.linspace(0.2, 0.4)
    mtop = 173.03
    xi = 0.
    tree_level_tadpoles = False

    lower = np.array([call_pt(l, 0.5 * mtop, xi, tree_level_tadpoles) for l in lambda_hs])
    central = np.array([call_pt(l, mtop, xi, tree_level_tadpoles) for l in lambda_hs])
    upper = np.array([call_pt(l, 2. * mtop, xi, tree_level_tadpoles) for l in lambda_hs])

    # TC and gamma

    for i in range(2):
        ax[i].plot(lambda_hs, lower[:, i], color="red", ls="--", label="$Q = 1/2 m_t$")
        ax[i].plot(lambda_hs, central[:, i], color="red", ls="-", label="$Q = m_t$")
        ax[i].plot(lambda_hs, upper[:, i], color="red", ls=":", label="$Q = 2 m_t$")

    # Higgs masses

    ax[2].plot(lambda_hs, central[:, 2], color="red", ls="--", label="Tree-level")
    ax[2].plot(lambda_hs, central[:, 3], color="red", ls="-", label="One-loop")

    compare = np.loadtxt("compare.txt")
    ax[0].plot(compare[:, 0], compare[:, 1], c="grey", ls="--")
    ax[0].plot(compare[:, 0], compare[:, 2], c="grey", ls=":")
    ax[1].plot(compare[:, 0], compare[:, 3], c="grey", ls="--")
    ax[1].plot(compare[:, 0], compare[:, 4], c="grey", ls=":")

    ax[0].legend(fontsize="small")
    ax[2].legend(fontsize="small")
    ax[0].set_ylim(0., 180.)
    ax[1].set_ylim(0., 4.)
    ax[0].set_ylabel("$T_C$ (GeV)")
    ax[1].set_ylabel("$\gamma$")
    ax[2].set_ylabel("$m_h$ (GeV)")
    for a in ax:
        a.set_xlabel("$\lambda_{hs}$")

    plt.tight_layout()
    name = "tree_level_tadpoles.png" if tree_level_tadpoles else "one_loop_tadpoles.png"
    plt.savefig(name)

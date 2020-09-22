#!/usr/bin/env python
"""
Compare h-bar expansion against https://arxiv.org/pdf/1808.01098.pdf fig. 2 
"""

from subprocess import check_output, run
import matplotlib.pyplot as plt
import numpy as np
import subprocess


def call_pt(lambda_hs, Q):
    r = subprocess.run(["../bin/h_bar_expansion", str(lambda_hs), str(Q)], stdout=subprocess.PIPE, encoding='utf8')

    for line in r.stdout.split("\n"):
        if line.startswith("TC = "):
            TC = float(line.split("=")[-1].strip())
        if line.startswith("gamma = "):
            gamma = float(line.split("=")[-1].strip())

    return TC, gamma

if __name__ == "__main__":

    lambda_hs = np.linspace(0.2, 0.4)
    mtop = 172.
    lower = np.array([call_pt(l, 0.5 * mtop) for l in lambda_hs])
    central = np.array([call_pt(l, mtop) for l in lambda_hs])
    upper = np.array([call_pt(l, 2. * mtop) for l in lambda_hs])

    fig, ax = plt.subplots(1, 2)
    for i, a in enumerate(ax):
        a.plot(lambda_hs, lower[:, i], c="red", ls="--", label="$Q = 1/2 m_t$")
        a.plot(lambda_hs, central[:, i], c="red", ls="-", label="$Q = m_t$")
        a.plot(lambda_hs, upper[:, i], c="red", ls=":", label="$Q = 2 m_t$")
        a.set_xlabel("$\lambda_{hs}$")

    ax[0].legend()
    ax[0].set_ylim(0., 180.)
    ax[1].set_ylim(0., 4.)
    ax[0].set_ylabel("$T_C$")
    ax[1].set_ylabel("$\gamma$")

    plt.savefig('h_bar_expansion.pdf')

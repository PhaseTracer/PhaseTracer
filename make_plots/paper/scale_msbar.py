#!/usr/bin/env python

from subprocess import run
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from call_pt import call_pt, transpose_list_dict
from style import style


if __name__ == "__main__":

    style()
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))

    lambda_hs = np.linspace(0.2, 0.4)
    mtop = 173.03
    xi = 0.
    tree_level_tadpoles = False
    tree_ewsb = False

    lower = transpose_list_dict([call_pt("msbar", l, 0.5 * mtop, xi, tree_level_tadpoles, tree_ewsb) for l in lambda_hs])
    central = transpose_list_dict([call_pt("msbar", l, mtop, xi, tree_level_tadpoles, tree_ewsb) for l in lambda_hs])
    upper = transpose_list_dict([call_pt("msbar", l, 2. * mtop, xi, tree_level_tadpoles, tree_ewsb) for l in lambda_hs])

    # TC and gamma  

    for i, n in enumerate(["TC", "gamma"]):
        ax[i].plot(lambda_hs, lower[n], color="red", ls="--", label="$Q = 1/2 m_t$")
        ax[i].plot(lambda_hs, central[n], color="red", ls="-", label="$Q = m_t$")
        ax[i].plot(lambda_hs, upper[n], color="red", ls=":", label="$Q = 2 m_t$")

    # ax[0].set_ylim(0., 180.)
    # ax[1].set_ylim(0., 4.)
    ax[0].set_ylabel("$T_C$ (GeV)")
    ax[1].set_ylabel("$\gamma$")

    for a in ax:
        a.set_xlabel("$\lambda_{hs}$")
        a.legend(fontsize="small")
    
    plt.suptitle(r"$\overline{\text{MS}}$")
    plt.tight_layout()
    plt.savefig("scale_msbar.pdf")

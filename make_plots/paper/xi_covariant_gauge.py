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
    tree_level_tadpoles = False
    tree_ewsb = False

    lower = transpose_list_dict([call_pt("covariant_gauge", l, mtop, 0., tree_level_tadpoles, tree_ewsb) for l in lambda_hs])
    central = transpose_list_dict([call_pt("covariant_gauge", l, mtop, 1., tree_level_tadpoles, tree_ewsb) for l in lambda_hs])
    upper = transpose_list_dict([call_pt("covariant_gauge", l, mtop, 2., tree_level_tadpoles, tree_ewsb) for l in lambda_hs])

    # TC and gamma  

    for i, n in enumerate(["TC", "gamma"]):
        ax[i].plot(lambda_hs, lower[n], color="red", ls="--", label=r"$\xi = 0$")
        ax[i].plot(lambda_hs, central[n], color="red", ls="-", label=r"$\xi = 1$")
        ax[i].plot(lambda_hs, upper[n], color="red", ls=":", label=r"$\xi = 2$")

    # ax[0].set_ylim(0., 180.)
    # ax[1].set_ylim(0., 4.)
    ax[0].set_ylabel("$T_C$ (GeV)")
    ax[1].set_ylabel("$\gamma$")

    for a in ax:
        a.set_xlabel("$\lambda_{hs}$")
        a.legend(fontsize="small")
    
    plt.suptitle("Covariant gauge")
    plt.tight_layout()
    plt.savefig("xi_covariant_gauge.pdf")

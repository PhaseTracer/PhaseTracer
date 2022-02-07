"""
Line plots for paper
====================
"""

import click
import matplotlib.pyplot as plt
import numpy as np

from plot_fun import fun_gamma, fun_diff
from style import style


def plot_for_1d(ax, data, x_num, label, column):
    sel = (data[:, 3] > 0) & (fun_gamma(data) > 0)
    x = data[:, x_num]
    TC = data[:, 4]
    gamma = fun_gamma(data)
    ax[0, column].plot(x[sel], TC[sel], label=label, alpha=1)
    ax[1, column].plot(x[sel], gamma[sel], label=label, alpha=1)

def ax_styling(ax):
    for a in ax.flatten():
        a.grid(axis='x', alpha=0.75)
        a.grid(axis='y', alpha=0.75)

    ax[0, 0].set_ylim(40, 170)
    ax[1, 0].set_ylim(0, 5)

    ax[0, 0].set_ylabel(r"$T_C$ (GeV)")
    ax[1, 0].set_ylabel(r"Strength, $\gamma_{\rm EW}$")

def schemes():
    fig, ax = plt.subplots(2, 3, figsize=(13, 6), sharey="row", sharex="col")
    ax_styling(ax)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    names = {"default": r"$\overline{\rm MS}$",
             "OSlike": r"OS-like",
             "HT": r"HT",
             "PRM_woFS_0L": r"PRM"}

    for k, v in names.items():
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/lambda_hs_{k}.txt"), 2, v, 0)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/lambda_s_{k}.txt"), 1, v, 1)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/m_s_{k}.txt"), 0, v, 2)

    ax[1, 0].set_xlabel(r"$\lambda_{hs}$")
    ax[1, 1].set_xlabel(r"$\lambda_{s}$")
    ax[1, 2].set_xlabel(r"$m_s$ (GeV)")

    ax[1, 0].set_xlim(0.1, 0.4)
    ax[1, 1].set_xlim(0.04, 0.2)
    ax[1, 2].set_xlim(40, 100)

    ax[0, 0].legend(loc=3)

    return fig, ax

def daisy():
    fig, ax = plt.subplots(2, 3, figsize=(13, 6), sharey="row", sharex="col")
    ax_styling(ax)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    names = ["lambda_hs", "lambda_s", "m_s"]

    for ii, name in enumerate(names):
        data0 = np.loadtxt(f"../1d_bks/{name}_noDaisy.txt")
        data1 = np.loadtxt(f"../1d_bks/{name}_default.txt")
        data2 = np.loadtxt(f"../1d_bks/{name}_Parwani.txt")
        show_data1 = fun_diff(data1, data0, data0, show_gamma=False, sort=False)
        show_data2 = fun_diff(data2, data0, data0, show_gamma=False, sort=False)
        show_data3 = fun_diff(data2, data1, data1, show_gamma=False, sort=False)

        a = ax[0, ii]
        a.plot(show_data1[:, 2-ii], show_data1[:, 4], label=r"${\rm AE}-{\rm ND}$", alpha=1)
        a.plot(show_data2[:, 2-ii], show_data2[:, 4], label=r"${\rm PW}-{\rm ND}$", alpha=1)
        a.plot(show_data3[:, 2-ii], show_data3[:, 4], label=r"${\rm PW}-{\rm AE}$", alpha=1)

        show_data1 = fun_diff(data1, data0, data0, show_gamma=True, sort=False)
        show_data2 = fun_diff(data2, data0, data0, show_gamma=True, sort=False)
        show_data3 = fun_diff(data2, data1, data1, show_gamma=True, sort=False)

        a = ax[1, ii]
        a.plot(show_data1[:, 2-ii], show_data1[:, 4],
               label=r"${\rm AE}-{\rm ND}$", alpha=1)
        a.plot(show_data2[:, 2-ii], show_data2[:, 4],
               label=r"${\rm PW}-{\rm ND}$", alpha=1)
        a.plot(show_data3[:, 2-ii], show_data3[:, 4],
               label=r"${\rm PW}-{\rm AE}$", alpha=1)


    ax[0, 0].set_ylabel(r"$\Delta T_C$ (GeV)")
    ax[1, 0].set_ylabel(r"$\Delta \gamma_{\rm EW}$")

    ax[1, 0].set_xlabel(r"$\lambda_{hs}$")
    ax[1, 1].set_xlabel(r"$\lambda_{s}$")
    ax[1, 2].set_xlabel(r"$m_s$ (GeV)")

    ax[1, 0].set_xlim(0.1, 0.4)
    ax[1, 1].set_xlim(0.04, 0.2)
    ax[1, 2].set_xlim(40, 100)

    ax[0, 0].legend(loc=3)

    ax[0, 0].autoscale()
    ax[1, 0].autoscale()

    return fig, ax

def scale_line():
    fig, ax = plt.subplots(2, 3, figsize=(13, 6), sharey="row", sharex="col")
    ax_styling(ax)

    ax[0, 0].set_ylim(0, 60)
    ax[1, 0].set_ylim(0, 50)
    
    names = {"lowT_noRGE_woFS_mt": r"$(Q=m_t)$",
             "lowT_noRGE_woFS_05mt": r"$(Q=\frac{1}{2}m_t)$",
             "lowT_noRGE_woFS_2mt": r"$(Q=2m_t)$",
}

    for k, v in names.items():
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/lambda_hs_{k}.txt"), 2, v, 0)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/lambda_s_{k}.txt"), 1, v, 1)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/m_s_{k}.txt"), 0, v, 2)

    ax[1, 0].legend(loc=2)

    ax[1, 0].set_xlabel(r"$\lambda_{hs}$")
    ax[1, 1].set_xlabel(r"$\lambda_{s}$")
    ax[1, 2].set_xlabel(r"$m_s$ (GeV)")

#    ax[1, 0].set_xlim(0.1, 0.4)
#    ax[1, 1].set_xlim(0.04, 0.2)
#    ax[1, 2].set_xlim(40, 100)

    return fig, ax

def xi(zoom_in=False):
    fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharey='row', sharex='col')
    ax_styling(ax)
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    add_name = '_zoom_in' if zoom_in else ''

    names = {"MSbar"+add_name: r"$\overline{\rm MS}$ + \texttt{1l\_self\_energy}",
             "MSbar_1L_EWSB"+add_name: r"$\overline{\rm MS}$ + \texttt{1l\_tadpole}",
             "MSbar_no"+add_name: r"$\overline{\rm MS}$ + \texttt{catastrophe}",
             "PRM"+add_name: r"PRM + \texttt{1l\_tadpole}",
             "PRM_0L"+add_name: r"PRM"}

    for k, v in names.items():
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/Rxi_{k}.txt"), 10, v, 0)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/covariant_{k}.txt"), 10, v, 1)

    data = np.loadtxt("../1d_bks/Rxi_HT.txt")
    TC = data[:, 4]
    gamma = fun_gamma(data)
    ax[0, 0].plot([0, 10], [TC[0], TC[-1]], label="HT", alpha=1)
    ax[0, 1].plot([0, 10], [TC[0], TC[-1]], label="HT", alpha=1)
    ax[1, 0].plot([0, 10], [gamma[0], gamma[-1]], label="HT", alpha=1)
    ax[1, 1].plot([0, 10], [gamma[0], gamma[-1]], label="HT", alpha=1)

    ax[0, 0].set_title(r"$R_\xi$ gauge")
    ax[0, 1].set_title(r"Covariant gauge")

    handles, labels = ax[1, 1].get_legend_handles_labels()
    plt.figlegend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(0.115, 0.75, 0.795, 0.2))
    plt.subplots_adjust(top=0.7)

    ax[1, 0].set_xlabel(r"$\xi$")
    ax[1, 1].set_xlabel(r"$\xi_W=\xi_B$")

    for a in ax[1, :]:
        a.set_xlim(0, 10)
    ax[0, 0].set_ylim(80, 120)
    ax[1, 0].set_ylim(1.5, 2.5)

    return fig, ax

def xi_zoom_in():
    fig, ax = xi(True)
    for a in ax[1, :]:
        a.set_xlim(0, 0.2)
    ax[0, 0].set_ylim(110, 114)
    ax[1, 0].set_ylim(1.75, 2.01)
    return fig, ax

@click.command()
@click.argument("plot_type")
def make_plot(plot_type):
    style()
    fig, ax = globals()[plot_type]()
    plt.savefig(f"1d_{plot_type}.pdf")

if __name__ == "__main__":
    make_plot()

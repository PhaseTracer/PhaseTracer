"""
Line plots for paper
====================
"""

import click
import matplotlib.pyplot as plt
import numpy as np

from plot_fun import fun_gamma, fun_diff
from style import style

def plot_for_1d_sel(ax, data, x_num, label, column):
    sel_TC = (data[:, 3] > 0) &(data[:, 8] > 1) 
    sel_gamma = (fun_gamma(data) > 0)
    sel = sel_TC &  sel_gamma
    x = data[:, x_num]
    TC = data[:, 4]
    gamma = fun_gamma(data)
    ax[0, column].plot(x[sel], TC[sel], label=label, alpha=1)
    ax[1, column].plot(x[sel], gamma[sel], label=label, alpha=1)

def plot_for_1d(ax, data, x_num, label, column):
    sel_TC = (data[:, 3] >= 0) 
    sel_gamma = (fun_gamma(data) > 0)
    x = data[:, x_num]
    TC = data[:, 4]
    gamma = fun_gamma(data)
    ax[0, column].plot(x[sel_TC], TC[sel_TC], label=label, alpha=1)
    ax[1, column].plot(x[sel_gamma], gamma[sel_gamma], label=label, alpha=1)

def ax_styling(ax):
    for a in ax.flatten():
        a.grid(axis='x', alpha=0.75)
        a.grid(axis='y', alpha=0.75)

    ax[0, 0].set_ylim(40, 170)
    ax[1, 0].set_ylim(0, 5)

    ax[0, 0].set_ylabel(r"$T_c$ (GeV)")
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
        show_data1 = fun_diff(data1, data0, data0, norm=True, show_gamma=False, sort=False)
        show_data2 = fun_diff(data2, data0, data0, norm=True, show_gamma=False, sort=False)
        show_data3 = fun_diff(data1, data2, data1, norm=True, show_gamma=False, sort=False)

        a = ax[0, ii]
        a.plot(show_data1[:, 2-ii], show_data1[:, 4], label=r"${\rm AE}-{\rm ND}$", alpha=1)
        a.plot(show_data2[:, 2-ii], show_data2[:, 4], label=r"${\rm PW}-{\rm ND}$", alpha=1)
        a.plot(show_data3[:, 2-ii], show_data3[:, 4], label=r"${\rm AE}-{\rm PW}$", alpha=1)

        show_data1 = fun_diff(data1, data0, data0, norm=True, show_gamma=True, sort=False)
        show_data2 = fun_diff(data2, data0, data0, norm=True, show_gamma=True, sort=False)
        show_data3 = fun_diff(data1, data2, data1, norm=True, show_gamma=True, sort=False)

        a = ax[1, ii]
        a.plot(show_data1[:, 2-ii], show_data1[:, 4],
               label=r"${\rm AE}-{\rm ND}$", alpha=1)
        a.plot(show_data2[:, 2-ii], show_data2[:, 4],
               label=r"${\rm PW}-{\rm ND}$", alpha=1)
        a.plot(show_data3[:, 2-ii], show_data3[:, 4],
               label=r"${\rm AE}-{\rm PW}$", alpha=1)


    ax[0, 0].set_ylabel(r"$\Delta_{\rm daisy} T_c/T_c$")
    ax[1, 0].set_ylabel(r"$\Delta_{\rm daisy} \gamma/\gamma$")

    ax[1, 0].set_xlabel(r"$\lambda_{hs}$")
    ax[1, 1].set_xlabel(r"$\lambda_{s}$")
    ax[1, 2].set_xlabel(r"$m_s$ (GeV)")

    ax[1, 0].set_xlim(0.1, 0.4)
    ax[1, 1].set_xlim(0.04, 0.2)
    ax[1, 2].set_xlim(40, 100)

    ax[1, 0].set_ylim(-0.05, 0.1)

    ax[1, 0].legend(loc=1)

    ax[0, 0].autoscale()
#    ax[1, 0].autoscale()

    return fig, ax

def scale_line():
    fig, ax = plt.subplots(2, 3, figsize=(13, 6), sharey="row", sharex="col")
    ax_styling(ax)

    ax[0, 0].set_ylim(0, 60)
    ax[1, 0].set_ylim(0, 30)
    
    names = {"lowT_05mt": r"$(Q=\frac{1}{2}m_t)$",
             "lowT_mt": r"$(Q=m_t)$",
             "lowT_2mt": r"$(Q=2m_t)$",
}

    for k, v in names.items():
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/lambda_hs_{k}.txt"), 2, v, 0)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/lambda_s_{k}.txt"), 1, v, 1)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/m_s_{k}.txt"), 0, v, 2)

    ax[1, 0].legend(loc=2)

    ax[1, 0].set_xlabel(r"$\lambda_{hs}$")
    ax[1, 1].set_xlabel(r"$\lambda_{s}$")
    ax[1, 2].set_xlabel(r"$m_s$ (GeV)")

    ax[1, 0].set_xlim(0.363, 0.38)
    ax[1, 1].set_xlim(0.045, 0.051)
    ax[1, 2].set_xlim(43, 46)

    return fig, ax

def xi(zoom_in=False):
    fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharey='row', sharex='col')
    ax_styling(ax)
    plt.subplots_adjust(hspace=0.3, wspace=0.1)

    add_name = '_zoom_in' if zoom_in else ''

    names = {"MSbar"+add_name: r"$\overline{\rm MS}$ + \texttt{GC\_SelfEnergy\_Sol}",
             "MSbar_1L_EWSB"+add_name: r"$\overline{\rm MS}$ + \texttt{GC\_Tadpole\_Sol}",
             "MSbar_no"+add_name: r"$\overline{\rm MS}$"}

    for k, v in names.items():
        plot_for_1d_sel(ax, np.loadtxt(f"../1d_bks/Rxi_{k}.txt"), 10, v, 0)
        plot_for_1d_sel(ax, np.loadtxt(f"../1d_bks/covariant_{k}.txt"), 10, v, 1)


    names = {"PRM"+add_name: r"PRM + \texttt{1LHiggs\_1LTad}",
             "PRM_0L"+add_name: r"PRM"}

    for k, v in names.items():
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/Rxi_{k}.txt"), 10, v, 0)
        plot_for_1d(ax, np.loadtxt(f"../1d_bks/covariant_{k}.txt"), 10, v, 1)

    data = np.loadtxt("../1d_bks/Rxi_HT.txt")
    TC = data[:, 4]
    gamma = fun_gamma(data)
    ax[0, 0].plot([0, 100], [TC[0], TC[-1]], label="HT", alpha=1)
    ax[0, 1].plot([0, 100], [TC[0], TC[-1]], label="HT", alpha=1)
    ax[1, 0].plot([0, 100], [gamma[0], gamma[-1]], label="HT", alpha=1)
    ax[1, 1].plot([0, 100], [gamma[0], gamma[-1]], label="HT", alpha=1)

    ax[0, 0].set_title(r"$R_\xi$ gauge")
    ax[0, 1].set_title(r"Covariant gauge")

    handles, labels = ax[1, 1].get_legend_handles_labels()

    if zoom_in:
        handles = handles[:-2]
        labels = labels[:-2]

    plt.figlegend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(0.115, 0.75, 0.795, 0.2))
    plt.subplots_adjust(top=0.7)

    ax[1, 0].set_xlabel(r"$\xi$")
    ax[1, 1].set_xlabel(r"$\xi$")

    for a in ax[1, :]:
        a.set_xlim(0, 60)
    ax[0, 0].set_ylim(30, 150)
    ax[1, 0].set_ylim(1.5, 5)

    return fig, ax

def xi_zoom_in():
    fig, ax = xi(True)
    for a in ax[1, :]:
        a.set_xlim(0, 1.)
    ax[0, 0].set_ylim(97.5, 117.5)
    ax[1, 0].set_ylim(1.70, 2.01)
    return fig, ax

@click.command()
@click.argument("plot_type")
def make_plot(plot_type):
    style()
    fig, ax = globals()[plot_type]()
    plt.savefig(f"1d_{plot_type}.pdf")

if __name__ == "__main__":
    make_plot()

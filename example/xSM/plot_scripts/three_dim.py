"""
3d scatter plot of change in TC and gamma and histograms
========================================================
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import click

from plot_fun import fun_gamma_line
from style import style


data_default = np.loadtxt(
    "../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_mt.txt")
data_mu_05 = np.loadtxt(
    "../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_05mt.txt")
data_mu_2 = np.loadtxt(
    "../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_2mt.txt")
data_xi_0 = np.loadtxt("../random_scan_results/xSM_MSbarxi0.txt")
data_xi_25 = np.loadtxt("../random_scan_results/xSM_MSbarxi25.txt")


max_num_cmap = matplotlib.cm.get_cmap('rainbow', 2)


def data():
    """
    @brief Parses data from disk
    """
    data_set = [data_default, data_mu_05, data_mu_2, data_xi_0, data_xi_25]
    len_data = len(data_default)
    for i, d in enumerate(data_set):
        if len(d) != len_data:
            raise RuntimeError("Length of data file " +
                               str(i) + " is wrong.")

    data_diff = []
    data_gamma = []
    for ii in range(len_data):
        ms = data_default[ii][0]
        lambda_s = data_default[ii][1]
        lambda_hs = data_default[ii][2]
        Tc_default = data_default[ii][4]
        gamma_default = fun_gamma_line(data_default[ii])
        flag_sel = True
        for jj in range(len(data_set)):
            if abs(data_set[jj][ii][0] - ms) > 0.01:
                print(data_set[jj][ii][0], ms)
                raise RuntimeError("Content of data file " +
                                   str(jj) + " is wrong.")

            vsjj = data_set[jj][ii][8]
            gammajj = fun_gamma_line(data_set[jj][ii])
            if gammajj <= 0 or vsjj < 10:
                flag_sel = False

            if flag_sel:
                d_xi = abs(data_xi_25[ii][4] - data_xi_0[ii][4])
                d_scale = abs(data_mu_05[ii][4] - data_mu_2[ii][4])
                gamma_xi = abs(fun_gamma_line(
                    data_xi_25[ii]) - fun_gamma_line(data_xi_0[ii]))
                gamma_scale = abs(fun_gamma_line(
                    data_mu_05[ii]) - fun_gamma_line(data_mu_2[ii]))

                d_set = [d_scale, d_xi]
                data_diff.append([ms, lambda_s, lambda_hs, np.where(
                    d_set == np.max(d_set))[0][0], max(d_set), Tc_default])

                gamma_set = [gamma_scale, gamma_xi]
                data_gamma.append(
                    [ms, lambda_s, lambda_hs, 0, max(gamma_set), gamma_default])

    data_diff.sort(key=(lambda x: -x[4]))
    diff = np.array(data_diff)
    diff_gamma = np.array(data_gamma)
    return diff, diff_gamma


def add_max_num_legend(ax, **kwargs):
    """
    @brief Add legend for source of maximum change
    """
    unique = [0, 1]
    colors = max_num_cmap(unique)
    patch = [matplotlib.patches.Rectangle(
        (0, 0), 1, 1, fc=c, ec=c) for c in colors]
    labels = ["Scale", "Gauge"]
    ax.legend(patch, labels, framealpha=0.95,
              title="Greatest impact on $T_c$", **kwargs)


def scatter_max_num(ax, diff, diff_gamma, nn):
    """
    @brief Scatter plot showing source of maximum change
    """
    ax.grid(axis='x', alpha=0.75)
    ax.grid(axis='y', alpha=0.75)
    style = dict(c=abs(diff[:, 3]), cmap=max_num_cmap, alpha=0.5,
                 edgecolor='none', rasterized=True, s=10)

    if nn == 0:
        ax.scatter(diff[:, 5], diff[:, 4], **style)
        ax.set_xlabel(r'$T_c$')
        ax.set_ylabel(r"$\max |\Delta T_c| $ (GeV)")
        add_max_num_legend(ax)

    elif nn == 1:
        ax.scatter(diff[:, 5], diff[:, 4] / diff[:, 5], **style)
        ax.set_xlabel(r'$T_c$')
        ax.set_ylabel(r'$\max |\Delta T_c| / T_c$')
    else:
        ax.scatter(diff_gamma[:, 5], diff_gamma[:, 4], **style)
        ax.set_xlabel(r'$\gamma_{\rm EW}$')
        ax.set_ylabel(r'$\max |\Delta \gamma_{\rm EW}|$')


def scatter(ax, diff, nx, ny, title, labels, cbar=False):
    """
    @brief Scatter physical quantity, e.g. TC or gamma
    """
    s = ax.scatter(diff[:, nx], diff[:, ny], c=diff[:, 4],
                   vmax=30, cmap="rainbow", edgecolor='none',
                   s=10, alpha=1, rasterized=True)
    ax.set_title(title)
    ax.set_xlabel(labels[nx])
    ax.set_ylabel(labels[ny])
    xmin = min(diff[:, nx])
    xmax = max(diff[:, nx])
    ymin = min(diff[:, ny])
    ymax = max(diff[:, ny])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    if cbar:
        fig = plt.gcf()
        fig.subplots_adjust(right=0.9, hspace=0.3, wspace=0.3, bottom=0.125)
        cbar_ax = fig.add_axes([0.915, 0.57, 0.02, 0.3])
        fig.colorbar(s, cax=cbar_ax)


def hist_plot(axs, nn, diff):
    """
    @brief Histogram of data
    """
    show_gamma = nn == 1
    ax = axs[nn, 0]
    ax.grid(axis='y', alpha=0.75)
    ax.hist(x=diff[:, 4], bins=30, color='seagreen',
            log=True, alpha=1, rwidth=0.85)
    if show_gamma:
        ax.set_xlabel(r'$\max |\Delta \gamma_{\rm EW}|$')
    else:
        ax.set_xlabel(r"$\max |\Delta T_c|$ (GeV)")

    ax.set_ylabel("Number of samples")

    ax = axs[nn, 1]
    ax.grid(axis='y', alpha=0.75)
    ax.hist(x=diff[:, 4]/diff[:, 5], bins=30,
            color='seagreen', log=True, alpha=1, rwidth=0.85)
    if show_gamma:
        ax.set_xlabel(r'$\Delta_{\rm max} |\gamma_{\rm EW}|/\gamma_{\rm EW}$')
    else:
        ax.set_xlabel(r"$\Delta_{\rm max}|T_c|/T_c $")
    ax.set_ylabel("Number of samples")


@click.command()
@click.argument("plot_type")
def make_plot(plot_type):

    style()
    diff, diff_gamma = data()

    if plot_type == "scatter":
        fig, axs = plt.subplots(2, 3, figsize=(15, 10))
        title = r"$\max |\Delta T_c| $ (GeV)"
        labels = [r'$M_s$ (GeV)', r'$\lambda_{S}$', r'$\lambda_{hs}$']

        scatter(axs[0, 0], diff, 1, 2, title, labels)
        scatter(axs[0, 1], diff, 0, 1, title, labels)
        scatter(axs[0, 2], diff, 0, 2, title, labels, cbar=True)
        scatter_max_num(axs[1, 0], diff, diff_gamma, 0)
        scatter_max_num(axs[1, 1], diff, diff_gamma, 1)
        scatter_max_num(axs[1, 2], diff, diff_gamma, 2)
        plt.savefig('3d_mu_scatter.pdf')
    else:
        fig, axs = plt.subplots(2, 2, figsize=(10, 10))
        hist_plot(axs, 0, diff)
        hist_plot(axs, 1, diff_gamma)
        fig.tight_layout()
        plt.savefig('3d_mu_hist.pdf')


if __name__ == "__main__":
    make_plot()

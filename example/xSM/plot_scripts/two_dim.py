"""
2d plots for paper
==================
"""

import click
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from plot_fun import fun_diff, fun_gamma_line
from three_dim import add_max_num_legend, max_num_cmap
from style import style


@click.command()
@click.argument("plot_type")
@click.argument("show")
def make_plot(plot_type, show):

    style()

    plot_xi = plot_type == "xi"
    plot_scale = plot_type == "scale"
    plot_scheme = plot_type == "scheme"
    plot_daisy = plot_type == "daisy"

    plot_max = plot_type == "max"

    show_deltagamma = show == "deltagamma"
    show_deltaT = show == "deltaT"
    show_max_num = show == "max_num"

    if plot_xi:
        figure_name = "xi"
    if plot_scale:
        figure_name = "scale"
    if plot_scheme:
        figure_name = "scheme"
    if plot_daisy:
        figure_name = "daisy"
    if plot_max:
        figure_name = "max"
        if show_max_num:
            figure_name = "max_num"
            show_deltaT = True

    if show_deltaT:
        if plot_xi:
            title = r"$T_c{(\xi=25)}-T_c{(\xi=0)}$ (GeV)"
        if plot_scale:
            title = r"$T_c{(Q=2m_t)}-T_c{(Q=\frac{1}{2}m_t)}$ (GeV)"
        if plot_scheme:
            title = r"$T_c(\overline{\rm MS})-T_c({\text{OS-like}})$ (GeV)"
        if plot_daisy:
            if show_deltagamma:
                title = r"$\left|\gamma_{\rm EW}^{\rm AE}-\gamma_{\rm EW}^{\rm PW}|/\gamma_{\rm EW}^{\rm AE}\right|$"
            else:
                title = r"$\left|T_c^{\rm PW}-T_c^{\rm AE}\right|$ (GeV)"
        if plot_max:
            title = r"$\max |\Delta T_c| $ (GeV)"
            if show_max_num:
                title = r""
        cm = 'rainbow'
    else:
        title = r"$T_c{(\xi=0)}$ (GeV)"
        cm = 'cool'

    labels = [r'$M_s$ (GeV)', r'$\lambda_{S}$', r'$\lambda_{hs}$']

    def make_plot(ax, par, cbar=False, legend=False):
        #############################
        if plot_xi:
            data2 = np.loadtxt("../2d_scan/"+par+"_default.txt")
            data1 = np.loadtxt("../2d_scan/"+par+"_xi25.txt")
            show_data = fun_diff(data1, data2, data1,
                                 show_gamma=show_deltagamma)
            print("max=",max(show_data[:, 4]))
            print("min=",min(show_data[:, 4]))
            vmax = 10
        if plot_scale:
            data2 = np.loadtxt("../2d_scan/"+par+"_05mt.txt")
            data1 = np.loadtxt("../2d_scan/"+par+"_2mt.txt")
            show_data = fun_diff(data1, data2, data1,
                                 show_gamma=show_deltagamma)
        if plot_scheme:
            data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
            data2 = np.loadtxt("../2d_scan/"+par+"_OSlike.txt")
            vmax = 6
            show_data = fun_diff(data1, data2, data1, show_gamma=show_deltagamma,
                                 norm=False, use_abs=False, gamma_min=0, sort=True)
        if plot_daisy:
            data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
            data2 = np.loadtxt("../2d_scan/"+par+"_Parwani.txt")
            vmax = 20
            show_data = fun_diff(data1, data2, data1,
                                 show_gamma=show_deltagamma)

        if plot_max:
            vmax = 16

            data_default = np.loadtxt("../2d_scan/"+par+"_default.txt")
            data_xi3 = np.loadtxt("../2d_scan/"+par+"_xi25.txt")
            data_05mt = np.loadtxt("../2d_scan/"+par+"_05mt.txt")
            data_2mt = np.loadtxt("../2d_scan/"+par+"_2mt.txt")

            data_set = [data_default, data_xi3, data_05mt, data_2mt]

            len_data = len(data_default)
            for ii in range(len(data_set)):
                if len(data_set[ii]) != len_data:
                    print("Length of data file " +
                          par + str(ii) + " is wrong.")
                    sys.exit()

            data_diff = []
            for ii in range(len_data):
                ms = data_default[ii][0]
                lambda_s = data_default[ii][1]
                lambda_hs = data_default[ii][2]
                flag_sel = True
                for jj in range(len(data_set)):
                    if abs(data_set[jj][ii][0] - ms) > 0.01:
                        print(data_set[jj][ii][0], ms)
                        print("Content of data file " +
                              par + str(jj) + " is wrong.")
                        sys.exit()

                    TCjj = data_set[jj][ii][4]
                    vsjj = data_set[jj][ii][8]
                    gammajj = fun_gamma_line(data_set[jj][ii])
                    if gammajj <= 0 or vsjj < 10:
                        flag_sel = False

                if flag_sel:
                    if show_deltaT:
                        d_xi = abs(data_default[ii][4] - data_xi3[ii][4])
                        d_scale = abs(data_05mt[ii][4] - data_2mt[ii][4])
                    else:
                        d_xi = abs(fun_gamma_line(
                            data_default[ii]) - fun_gamma_line(data_xi3[ii]))
                        d_scale = abs(fun_gamma_line(
                            data_05mt[ii]) - fun_gamma_line(data_2mt[ii]))

                    d_set = [d_scale, d_xi]
                    data_diff.append([ms, lambda_s, lambda_hs, np.where(
                        d_set == np.max(d_set))[0][0], max(d_set)])
            show_data = np.array(data_diff)

        # 0=lambda_hs, 1=lambda_s, 2=ms
        if par == "lhs_ls":
            nx = 1
            ny = 2
            label = r"$M_s=65$ GeV"

        if par == "ms_lhs":
            nx = 0
            ny = 2
            label = r"$\lambda_{S}=0.1$"

        if par == "ms_ls":
            nx = 0
            ny = 1
            label = r"$\lambda_{hs}=0.3$"

        xmin = min(show_data[:, nx])
        xmax = max(show_data[:, nx])
        ymin = min(show_data[:, ny])
        ymax = max(show_data[:, ny])

        if show_deltaT:
            if plot_daisy:
                if show_deltagamma:
                    map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=abs(
                        show_data[:, 4]), cmap=cm, edgecolor='none', s=5, vmin=0, vmax=.2, alpha=1, rasterized=True)
                else:
                    map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=show_data[:, 4],
                                      cmap=cm, edgecolor='none', s=5, vmin=-0.1, vmax=4, alpha=1, rasterized=True)
            elif plot_scheme:
                map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=show_data[:, 4],
                                  cmap=cm, edgecolor='none', s=5, vmin=-1, vmax=10, alpha=1, rasterized=True)
            elif plot_max and show_max_num:
                map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=abs(
                    show_data[:, 3]), cmap=max_num_cmap, edgecolor='none', s=5, alpha=1, rasterized=True)
            elif plot_max and not show_max_num:
                map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=show_data[:, 4],
                                  cmap=cm, edgecolor='none', s=5, vmin=5, vmax=30, alpha=1, rasterized=True)
            elif plot_scale:
                map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=show_data[:, 4],
                                  cmap=cm, edgecolor='none', s=5, vmin=8, vmax=16, alpha=1, rasterized=True)
                print("max",max(show_data[:, 4]))
                print("par",show_data[np.where(show_data[:, 4] == np.max(show_data[:, 4]))[0][0]])
            elif plot_xi:
                map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=show_data[:, 4],
                                  cmap=cm, edgecolor='none', s=5, vmin=-73, vmax=55, alpha=1, rasterized=True)
                if par == "ms_ls":
                  xmin = 10.0
                  xmax = 90.7035
                  ymin = 0.01
                  ymax = 0.3
                if par == "ms_lhs":
                  xmin = 10.0
                  xmax = 115.025
                  ymin = 0.1
                  ymax = 0.5
        else:
            data2 = np.loadtxt("../2d_scan/"+par+"_default.txt")
            data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
            show_data = fun_diff(data1, data2, data1,
                                 show_gamma=show_deltagamma)
            xmin = min(show_data[:, nx])
            xmax = max(show_data[:, nx])
            ymin = min(show_data[:, ny])
            ymax = max(show_data[:, ny])
            map1 = ax.scatter(show_data[:, nx], show_data[:, ny], c=show_data[:, 5],
                              cmap=cm, s=2, vmax=150, alpha=1, rasterized=True)
            print("xmin",xmin,"xmax",xmax)
            print("ymin",ymin,"ymax",ymax)

        ax.set_xlabel(labels[nx])
        ax.set_ylabel(labels[ny])
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_title(title, fontsize=14)
        ax.text(xmin+0.1*(xmax-xmin), ymax-0.1*(ymax-ymin), label)

        if plot_max and show_max_num and legend:
            add_max_num_legend(ax, loc="lower right")

        if cbar:
            fig = plt.gcf()
            fig.subplots_adjust(right=0.9, wspace=0.3, bottom=0.125)

            if not show_max_num:
                cbar_ax = fig.add_axes([0.915, 0.15, 0.02, 0.7])
                fig.colorbar(map1, cax=cbar_ax)
                if plot_daisy and show_deltagamma:
                    clb.set_ticks([0, .05, .1, .15, .2])
                    clb.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    make_plot(axs[0], par="lhs_ls", legend=True)
    make_plot(axs[1], par="ms_ls")
    make_plot(axs[2], par="ms_lhs", cbar=True)

    figname = '2d_scan_'+figure_name

    if show_deltaT:
        if plot_daisy and show_deltagamma:
            figname += '_deltagamma'
        else:
            figname += '_deltaT'
    else:
        figname += '_T'

    plt.savefig(figname+'.pdf')


if __name__ == "__main__":
    make_plot()

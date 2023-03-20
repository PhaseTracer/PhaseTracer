"""
Band plots for paper
====================

"plot_scale" in the first place need to be manually modified.

"""

import click
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from plot_fun import fun_gamma
from style import style


style()


plot_scale = True

def selection(for_TC, x_num, data, PRM=False):
    x = data[:, x_num]

    if for_TC:
      if plot_scale:
        sel = ( (data[:,3]>0) | (PRM&(data[:, 3]>=0)) ) & (data[:, 4] < 200.) 
        y = data[:, 4]
      else:
        sel = ( (data[:,3]>0) | (PRM&(data[:, 3]>=0)) ) & (data[:, 4] < 200.) 
        Q_sel = ( abs(data[:,9]/data[:,4] - 3.1415926*2) < 0.01 ) | ( abs(data[:,9]/data[:,4] - 0.5) < 0.01 ) | ( abs(data[:,9]/data[:,4] - 1) < 0.01 ) 
        sel = sel & Q_sel
        y = data[:, 4]
    else:
        sel = ( (data[:,3]>0) | (PRM &(data[:, 3]>=0)) )
        y = fun_gamma(data)
        # allow points with one minima in PRM method
        if PRM:
            where = y == 0
            y[where] = np.maximum(abs(data[:, 7][where]), abs(data[:, 5][where])) / data[:, 4][where]
    return x[sel], y[sel]

def line_for_1d(axs, data, x_num, color, column, label="_nolegend_", linestyle="-", alpha=1, PRM=False):
    x_TC, TC = selection(True, x_num, data, PRM)
    x_gamma, gamma = selection(False, x_num, data, PRM)

    ax = axs[0, column]
    ax.plot(x_TC, TC, color=color, linestyle=linestyle, alpha=alpha, label=label)

    if plot_scale:
      ax = axs[2, column]
      ax.plot(x_gamma, gamma, color=color, linestyle=linestyle, alpha=alpha, label=label)

    fTC = interp1d(x_TC, TC, kind='linear', bounds_error=False, fill_value='extrapolate')
    fgamma = interp1d(x_gamma, gamma, kind='linear', bounds_error=False, fill_value='extrapolate')
    return fTC, fgamma

def interpolate(x1, y1, x2, y2):
    if min(x1) < min(x2):
        x2 = np.append(x1[0], x2)
        y2 = np.append(y1[0], y2)
    elif min(x1) > min(x2):
        x1 = np.append(x2[0], x1)
        y1 = np.append(y2[0], y1)
    if max(x1) > max(x2):
        x2 = np.append(x2, x1[-1])
        y2 = np.append(y2, y1[-1])
    elif max(x1) < max(x2):
        x1 = np.append(x1, x2[-1])
        y1 = np.append(y1, y2[-1])

    f1 = interp1d(x1, y1, kind='linear')
    f2 = interp1d(x2, y2, kind='linear')
    return f1, f2

def lucd(x, f1, f2, f):
    f1x = f1(x)
    f2x = f2(x)
    fx = f(x)
    for i in range(len(fx)):
      if f2x[i] > f1x[i]:
        if f2x[i] < fx[i]:
          f2x[i] = fx[i]
        if f1x[i] > fx[i]:
          f1x[i] = fx[i]
      else:
        if f2x[i] > fx[i]:
          f2x[i] = fx[i]
        if f1x[i] < fx[i]:
          f1x[i] = fx[i]
          
    return f1x, f2x, fx, f1x - f2x

def range_for_1d(axs, data1, data2, x_num, label, column, color, fTC, fgamma, PRM=False):

    x_TC1, TC1 = selection(True, x_num, data1, PRM)
    x_TC2, TC2 = selection(True, x_num, data2, PRM)

    x_gamma1, gamma1 = selection(False, x_num, data1, PRM)
    x_gamma2, gamma2 = selection(False, x_num, data2, PRM)

    fTC1, fTC2 = interpolate(x_TC1, TC1, x_TC2, TC2)
    fgamma1, fgamma2 = interpolate(x_gamma1, gamma1, x_gamma2, gamma2)

    x_TC = np.linspace(min(min(x_TC1), min(x_TC2)), max(max(x_TC1),max(x_TC2)), num=100, endpoint=True)
    x_gamma = np.linspace(min(min(x_gamma1), min(x_gamma2)), max(max(x_gamma1), max(x_gamma2)), num=100, endpoint=True)

    lTC, uTC, cTC, dTC = lucd(x_TC, fTC1, fTC2, fTC)
    ax = axs[0,column]
    ax.fill_between(x_TC, lTC, uTC, color=color, alpha=0.3, linewidth=0, label=label)

    if plot_scale:
      lgamma, ugamma, cgamma, dgamma = lucd(x_gamma, fgamma1, fgamma2, fgamma)
      ax = axs[2,column]
      ax.fill_between(x_gamma, lgamma, ugamma, color=color, alpha=0.3, linewidth=0, label=label)

    x_TC = np.linspace(max(min(x_TC1), min(x_TC2)), min(max(x_TC1), max(x_TC2)), num=100, endpoint=True)
    lTC, uTC, cTC, dTC = lucd(x_TC, fTC1, fTC2, fTC)
    ax = axs[1,column]
    ax.plot(x_TC, dTC/cTC, color=color, label=label)

    if plot_scale:
      x_gamma = np.linspace(max(min(x_gamma1), min(x_gamma2)), min(max(x_gamma1), max(x_gamma2)), num=100, endpoint=True)
      lgamma, ugamma, cgamma, dgamma = lucd(x_gamma, fgamma1, fgamma2, fgamma)
      ax = axs[3,column]
      ax.plot(x_gamma, dgamma / cgamma, color=color, label=label)

#    ax = axs[0,column]
#    ax.fill_between(x_TC,fTC1(x_TC),fTC2(x_TC), color=color, alpha=0.3, linewidth=0, label=label)
#
#    ax = axs[2,column]
#    ax.fill_between(x_gamma,fgamma1(x_gamma),fgamma2(x_gamma), color=color, alpha=0.3, linewidth=0, label=label)
#
#    x_TC = np.linspace(max(min(x_TC1), min(x_TC2)), min(max(x_TC1), max(x_TC2)), num=100, endpoint=True)
#    ax = axs[1,column]
#    ax.plot(x_TC,(fTC1(x_TC)-fTC2(x_TC))/fTC(x_TC), color=color, label=label)
#
#    x_gamma = np.linspace(max(min(x_gamma1), min(x_gamma2)), min(max(x_gamma1), max(x_gamma2)), num=100, endpoint=True)
#    ax = axs[3,column]
#    ax.plot(x_gamma,(fgamma1(x_gamma)-fgamma2(x_gamma))/fgamma(x_gamma), color=color, label=label)



@click.command()
@click.argument("plot_type")
def make_plot(plot_type):

    if plot_type == "scale":
        plot_scale = True
    elif plot_type == "scaleT":
        plot_scale = False
    else:
        raise RuntimeError("unknown plot type")

    if plot_scale:
      fig, axs = plt.subplots(4, 3, figsize=(11, 12), sharey="row", sharex="col")
    else:
      fig, axs = plt.subplots(2, 3, figsize=(11, 6), sharey="row", sharex="col")
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

    for name in [["m_s", 0, 2], ["lambda_s", 1, 1], ["lambda_hs", 2, 0]]:

      if plot_scale:

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_mt.txt"), name[1], "purple", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_05mt.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_2mt.txt"),
                     name[1], r"$\overline{\rm MS}$ w AE", name[2], 'purple', fTC, fgamma)
                     
        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PW_mt.txt"), name[1], "g", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PW_05mt.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_PW_2mt.txt"),
                     name[1], r"$\overline{\rm MS}$ w PW", name[2], 'g', fTC, fgamma)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noD_mt.txt"), name[1], "b", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noD_05mt.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_noD_2mt.txt"),
                     name[1], r"$\overline{\rm MS}$ w/o daisy", name[2], 'b', fTC, fgamma)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noRGE_woFS_mt.txt"), name[1], "gray", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noRGE_woFS_05mt.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_noRGE_woFS_2mt.txt"),
                     name[1], r"$\overline{\rm MS}$ w/o RGE", name[2], 'gray', fTC, fgamma)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_mt.txt"), name[1], "r", name[2], PRM=True)
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_05mt.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_2mt.txt"),
                     name[1], r"PRM", name[2], 'r', fTC, fgamma, PRM=True)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_mt.txt"), name[1], "orange", name[2], PRM=True)
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_05mt.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_2mt.txt"),
                     name[1], r"PRM w/o RGE", name[2], 'orange', fTC, fgamma, PRM=True)
#        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_2mt.txt"), name[1], "purple", name[2], linestyle="--", alpha=0.4, PRM=True)

      else:
      
        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_T.txt"), name[1], "purple", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_05T.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_2piT.txt"),
                     name[1], r"$\overline{\rm MS}$ w AE", name[2], 'purple', fTC, fgamma)
                     
        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PW_T.txt"), name[1], "g", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PW_05T.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_PW_2piT.txt"),
                     name[1], r"$\overline{\rm MS}$ w PW", name[2], 'g', fTC, fgamma)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noD_T.txt"), name[1], "b", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noD_05T.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_noD_2piT.txt"),
                     name[1], r"$\overline{\rm MS}$ w/o daisy", name[2], 'b', fTC, fgamma)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noRGE_woFS_T.txt"), name[1], "gray", name[2])
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_noRGE_woFS_05T.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_noRGE_woFS_2piT.txt"),
                     name[1], r"$\overline{\rm MS}$ w/o RGE", name[2], 'gray', fTC, fgamma)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_T.txt"), name[1], "r", name[2], PRM=True)
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_05T.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_2piT.txt"),
                     name[1], r"PRM", name[2], 'r', fTC, fgamma, PRM=True)

        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_T.txt"), name[1], "orange", name[2], PRM=True)
        range_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_05T.txt"),
                     np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_2piT.txt"),
                     name[1], r"PRM w/o RGE", name[2], 'orange', fTC, fgamma, PRM=True)
#        fTC, fgamma = line_for_1d(axs, np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_2piT.txt"), name[1], "purple", name[2], linestyle="--", alpha=0.4, PRM=True)


    for ii in range(4 if plot_scale else 2):
      for jj in range(3):
        axs[ii, jj].grid(axis='x', alpha=0.5)
        axs[ii, jj].grid(axis='y', alpha=0.5)

        if ii == 0:
          axs[ii, jj].set_ylim(20, 170)
          axs[ii, 0].set_ylabel(r"$T_c$ (GeV)")
        elif ii == 2:
          axs[ii, jj].set_ylim(0, 6)
          axs[ii, 0].set_ylabel(r"$\gamma_{\rm EW}$")
        elif ii == 1:
          axs[ii, jj].set_ylim(-0.5, 0.5)
          axs[ii, 0].set_ylabel(r"$\Delta_Q T_c/ T_c$")
        elif ii == 3:
          axs[ii, jj].set_ylim(-1.5, 1)
          axs[ii, 0].set_ylabel(r"$\Delta_Q \gamma / \gamma$")

        if jj == 0:
          axs[-1, jj].set_xlabel(r"$\lambda_{hs}$")
          axs[ii, jj].set_xlim(0.1, 0.4)
        elif jj == 1:
          axs[-1, jj].set_xlabel(r"$\lambda_{s}$")
          axs[ii, jj].set_xlim(0.05, 0.2)
        else:
          axs[-1, jj].set_xlabel(r"$m_{s}$ (GeV)")
          axs[ii, jj].set_xlim(40, 100)


    if plot_scale:
      handles, labels = axs[1, 1].get_legend_handles_labels()
      plt.figlegend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.115, 0.75, 0.795, 0.2))
      plt.subplots_adjust(top=0.87)
    else:
      handles, labels = axs[1, 1].get_legend_handles_labels()
      plt.figlegend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.115, 0.80, 0.795, 0.2))
      plt.subplots_adjust(top=0.85)

    if plot_scale:
      plt.savefig('1d_scale.pdf')
    else:
      plt.savefig('1d_scale_T.pdf')

if __name__ == "__main__":
    make_plot()

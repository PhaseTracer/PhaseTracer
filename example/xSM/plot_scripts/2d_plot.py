"""
2d plots for paper
==================
"""

import click
import matplotlib.pyplot as plt
import numpy as np
from plot_fun import fun_diff, fun_gamma_line
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
        title=r"$\left|T_C{(\xi=3)}-T_C{(\xi=0)}\right|$ (GeV)"
      if plot_scale:
        title=r"$\left|T_C{(Q=\frac{1}{2}m_t)}-T_C{(Q=2m_t)}\right|$ (GeV)"
      if plot_scheme:
        title=r"$\left|T_C^{\overline{\rm MS}}-T_C^{\rm OS-like}\right|$ (GeV)"
      if plot_daisy:
        if show_deltagamma:
          title=r"$\left|\gamma_{\rm EW}^{\rm AE}-\gamma_{\rm EW}^{\rm PW}|/\gamma_{\rm EW}^{\rm AE}\right|$"
        else: 
          title=r"$\left|T_C^{\rm PW}-T_C^{\rm AE}\right|$ (GeV)"
      if plot_max:
        title=r"$\Delta T_C ^{\rm max}$ (GeV)"
        if show_max_num:
          title=r"Source of $\Delta T_C ^{\rm max}$ (GeV)"
      cm = 'rainbow'
    else:
      title=r"$T_C{(\xi=0)}$ (GeV)"
      cm = 'cool'

    labels=[r'$M_s$ (GeV)', r'$\lambda_{S}$', r'$\lambda_{hs}$']

    def make_plot(ax, par, cbar=False):
      #############################
      if plot_xi:
        data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
        data2 = np.loadtxt("../2d_scan/"+par+"_xi3.txt")
        vmax = 10
      if plot_scale:
        data1 = np.loadtxt("../2d_scan/"+par+"_05mt.txt")
        data2 = np.loadtxt("../2d_scan/"+par+"_2mt.txt")
        vmax = 16
      if plot_scheme:
        data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
        data2 = np.loadtxt("../2d_scan/"+par+"_OSlike.txt")
        vmax = 6
      if plot_daisy:
        data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
        data2 = np.loadtxt("../2d_scan/"+par+"_Parwani.txt")
        vmax = 20
      
      if plot_max:
        vmax = 16

        data_default = np.loadtxt("../2d_scan/"+par+"_default.txt")
        data_xi3 = np.loadtxt("../2d_scan/"+par+"_xi3.txt")
        data_05mt = np.loadtxt("../2d_scan/"+par+"_05mt.txt")
        data_2mt = np.loadtxt("../2d_scan/"+par+"_2mt.txt")
        data_Parwani = np.loadtxt("../2d_scan/"+par+"_Parwani.txt")
        data_OSlike = np.loadtxt("../2d_scan/"+par+"_OSlike.txt")
        
        data_set = [data_default, data_xi3, data_05mt, data_2mt, data_Parwani, data_OSlike]
        
        len_data = len(data_default)
        for ii in range(len(data_set)): 
          if len(data_set[ii]) != len_data:
            print("Length of data file is wrong.")
            sys.exit()
         
        data_diff=[]
        for ii in range(len_data):
          ms = data_default[ii][0]
          lambda_s = data_default[ii][1]
          lambda_hs = data_default[ii][2]
          flag_sel = True
          for jj in range(len(data_set)): 
            if data_set[jj][ii][0] != ms:
              print("Content of data file is wrong.")
              sys.exit()
            
            TCjj = data_set[jj][ii][4]
            vsjj = data_set[jj][ii][8]
            gammajj = fun_gamma_line(data_set[jj][ii])
            if gammajj<=0 or  vsjj<10:
              flag_sel = False
          
          if flag_sel:
            if show_deltaT:
              d_xi = abs(data_default[ii][4] - data_xi3[ii][4])
              d_scale = abs(data_05mt[ii][4] - data_2mt[ii][4])
              d_scheme = abs(data_default[ii][4] - data_OSlike[ii][4])
              d_daisy = abs(data_default[ii][4] - data_Parwani[ii][4])
            else: 
              d_xi = abs(fun_gamma_line(data_default[ii]) - fun_gamma_line(data_xi3[ii]))
              d_scale = abs(fun_gamma_line(data_05mt[ii]) - fun_gamma_line(data_2mt[ii]))
              d_scheme = abs(fun_gamma_line(data_default[ii]) - fun_gamma_line(data_OSlike[ii]))
              d_daisy = abs(fun_gamma_line(data_default[ii]) - fun_gamma_line(data_Parwani[ii]))
      
            d_set = [d_xi, d_scheme, d_daisy]
            data_diff.append([ms, lambda_s, lambda_hs, np.where(d_set==np.max(d_set))[0][0], max(d_set)])
        show_data = np.array(data_diff)
       
       
      
      # 0=lambda_hs, 1=lambda_s, 2=ms
      if par=="lhs_ls":    
        nx = 1
        ny = 2
        label = r"$M_s=65$ GeV"
        
      if par=="ms_lhs":
        nx = 0
        ny = 2
        label = r"$\lambda_{S}=0.1$"

      if par=="ms_ls":
        nx = 0
        ny = 1
        label = r"$\lambda_{hs}=0.3$"

      xmin = min(show_data[:,nx])
      xmax = max(show_data[:,nx])
      ymin = min(show_data[:,ny])
      ymax = max(show_data[:,ny])
        
      if show_deltaT:
        if plot_daisy:
          if show_deltagamma:
            map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=abs(show_data[:,4]), cmap=cm, edgecolor='none', s=5, vmin=0, vmax=.2, alpha=1, rasterized=True)
          else:
            map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=show_data[:,4], cmap=cm, edgecolor='none', s=5, vmin=-0.1, vmax=2, alpha=1, rasterized=True)
        elif plot_scheme:
          map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=-show_data[:,4], cmap=cm, edgecolor='none', s=5, vmin=-1, vmax=10, alpha=1, rasterized=True)
        elif plot_max and show_max_num:
          map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=abs(show_data[:,3]), cmap=cm, edgecolor='none', s=5, vmin=0, vmax=4, alpha=1, rasterized=True)
        else:
          map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=abs(show_data[:,4]), cmap=cm, edgecolor='none', s=5, vmax=vmax, alpha=1, rasterized=True)
      else:
        map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=show_data[:,5], cmap=cm, s=2, vmax=150, alpha=1, rasterized=True)
        
      ax.set_xlabel(labels[nx])
      ax.set_ylabel(labels[ny])
      ax.set_xlim(xmin,xmax)
      ax.set_ylim(ymin,ymax)
      ax.set_title(title, fontsize=14)
      ax.text(xmin+0.1*(xmax-xmin),ymax-0.1*(ymax-ymin),label)

      if cbar:
          fig = plt.gcf()
          fig.subplots_adjust(right=0.9, wspace=0.3, bottom=0.125)
          cbar_ax = fig.add_axes([0.915, 0.15, 0.02, 0.7])
          fig.colorbar(map1, cax=cbar_ax)
          if plot_daisy and show_deltagamma:
            clb.set_ticks([0, .05, .1, .15, .2])
            clb.set_ticklabels(['0%', '5%', '10%', '15%', '20%'])


    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    make_plot(axs[0], par="lhs_ls")
    make_plot(axs[1], par="ms_ls")
    make_plot(axs[2], par="ms_lhs", cbar=True)

    #fig.tight_layout()

    #plt.subplots_adjust(hspace=0.2)

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

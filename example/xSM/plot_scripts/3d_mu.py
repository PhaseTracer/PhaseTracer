import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_fun import fun_gamma, fun_diff, loaddata, fun_gamma_line
from style import style
import click

data_default = np.loadtxt("../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_mt.txt")
data_mu_05 = np.loadtxt("../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_05mt.txt")
data_mu_2 = np.loadtxt("../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_2mt.txt")
data_xi_0 = np.loadtxt("../random_scan_results/xSM_MSbarxi0.txt")
data_xi_25 = np.loadtxt("../random_scan_results/xSM_MSbarxi25.txt")
  
  
show_gamma = False
use_log = False

@click.command()
@click.argument("plot_type")
def make_plot(plot_type):

  style()
    
  data_set = [data_default, data_mu_05, data_mu_2, data_xi_0, data_xi_25]    
  len_data = len(data_default)
  for ii in range(len(data_set)): 
      if len(data_set[ii]) != len_data:
        print("Length of data file " +par + str(ii) +" is wrong.")
        sys.exit()
         
  data_diff=[]
  data_gamma=[]
  for ii in range(len_data):
      ms = data_default[ii][0]
      lambda_s = data_default[ii][1]
      lambda_hs = data_default[ii][2]
      Tc_default = data_default[ii][4]
      gamma_default = fun_gamma_line(data_default[ii])
      flag_sel = True
      for jj in range(len(data_set)): 
        if abs( data_set[jj][ii][0] - ms) > 0.01:
          print(data_set[jj][ii][0], ms)
          print("Content of data file " +par  + str(jj) +" is wrong.")
          sys.exit()
            
        TCjj = data_set[jj][ii][4]
        vsjj = data_set[jj][ii][8]
        gammajj = fun_gamma_line(data_set[jj][ii])
        if gammajj<=0 or  vsjj<10:
          flag_sel = False
          
        if flag_sel:
          d_xi = abs(data_xi_25[ii][4] - data_xi_0[ii][4])
          d_scale = abs(data_mu_05[ii][4] - data_mu_2[ii][4])
          gamma_xi = abs(fun_gamma_line(data_xi_25[ii]) - fun_gamma_line(data_xi_0[ii]))
          gamma_scale = abs(fun_gamma_line(data_mu_05[ii]) - fun_gamma_line(data_mu_2[ii]))      
      
          d_set = [d_scale, d_xi]
          data_diff.append([ms, lambda_s, lambda_hs, np.where(d_set==np.max(d_set))[0][0], max(d_set), Tc_default])

          gamma_set = [gamma_scale, gamma_xi]
          data_gamma.append([ms, lambda_s, lambda_hs, 0, max(gamma_set), gamma_default])

  data_diff.sort(key=(lambda x:-x[4])) 
  diff = np.array(data_diff)
  diff_gamma = np.array(data_gamma)
    
  if plot_type == "scatter":
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    vmax = 20
    size = 10
    title =  r"$\Delta_{\rm max}|T_c| $ (GeV)"
    labels=[r'$M_s$ (GeV)', r'$\lambda_{S}$', r'$\lambda_{hs}$']
      
      
    def make_plot(ax, nx, ny, cbar=False):
      map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,4], vmax = 30, cmap="rainbow", edgecolor='none', s=size, alpha=1, rasterized=True)
      ax.set_title(title)
      ax.set_xlabel(labels[nx])
      ax.set_ylabel(labels[ny])
      xmin = min(diff[:,nx])
      xmax = max(diff[:,nx])
      ymin = min(diff[:,ny])
      ymax = max(diff[:,ny])
      ax.set_xlim(xmin,xmax)
      ax.set_ylim(ymin,ymax)
    
      if cbar:
        fig = plt.gcf()
        fig.subplots_adjust(right=0.9, wspace=0.3, bottom=0.125)
            
        cbar_ax = fig.add_axes([0.915, 0.50, 0.02, 0.4])
        fig.colorbar(map, cax=cbar_ax)

    def no_color(ax, nn):
      ax.grid(axis='x', alpha=0.75)
      ax.grid(axis='y', alpha=0.75)
      if nn == 0 :
        ax.scatter(diff[:,5], diff[:,4], c='seagreen', alpha=1, edgecolor='none', rasterized=True)
        ax.set_xlabel(r'$T_c(Q=m_t, \xi=1)$')
        ax.set_ylabel(r"$\Delta_{\rm max}|T_c| $ (GeV)")
      elif nn == 1 :
        ax.scatter(diff[:,5], diff[:,4]/diff[:,5], c='seagreen', alpha=1, edgecolor='none', rasterized=True)
        ax.set_xlabel(r'$T_c(Q=m_t, \xi=1)$')
        ax.set_ylabel(r'$\Delta_{\rm max}|T_c|/T_c(Q=m_t, \xi=1)$')
      else:
        ax.scatter(diff_gamma[:,5], diff_gamma[:,4], c='seagreen', alpha=1, edgecolor='none', rasterized=True)
        ax.set_xlabel(r'$\gamma_{\rm EW}$')
        ax.set_ylabel(r'$\Delta_{\rm max} |\gamma_{\rm EW}|$')
#        ax.set_ylim(-5,0.1)
#        ax.set_xlim(-0.1,7)
            
    make_plot(axs[0,0], 1, 2)
    make_plot(axs[0,1], 0, 1)
    make_plot(axs[0,2], 0, 2, cbar=True)
    no_color(axs[1,0], 0)
    no_color(axs[1,1], 1)
    no_color(axs[1,2], 2)
    plt.savefig('3d_mu_scatter.pdf')
  
  else:
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    
    def hist_plot(nn, diff):
      show_gamma = (nn==1)
      ax = axs[nn,0]
      ax.grid(axis='y', alpha=0.75)
      ax.hist(x=diff[:,4], bins=30, color='seagreen', log=True, alpha=1, rwidth=0.85)
      if show_gamma:
        ax.set_xlabel(r'$\Delta_{\rm max} |\gamma_{\rm EW}|$')
      else:
        ax.set_xlabel(r"$\Delta_{\rm max}|T_c| $ (GeV)")
      
      ax.set_ylabel("Number of samples")
    
      ax = axs[nn,1]
      ax.grid(axis='y', alpha=0.75)
      ax.hist(x=diff[:,4]/diff[:,5], bins=30, color='seagreen', log=True, alpha=1, rwidth=0.85)
      if show_gamma:
        ax.set_xlabel(r'$\Delta_{\rm max} |\gamma_{\rm EW}|/\gamma_{\rm EW}$')
      else:
        ax.set_xlabel(r"$\Delta_{\rm max}|T_c|/T_c $")
      ax.set_ylabel("Number of samples")
    
    hist_plot(0,diff)
    hist_plot(1,diff_gamma)
    fig.tight_layout()
    plt.savefig('3d_mu_hist.pdf')
  
  plt.show()
 
if __name__ == "__main__":
    make_plot()
    

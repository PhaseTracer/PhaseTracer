import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_fun import fun_gamma, fun_diff, loaddata
from style import style
import click

data_mu_05 = np.loadtxt("../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_05mt.txt")
data_mu_1 = np.loadtxt("../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_mt.txt")
data_mu_2 = np.loadtxt("../random_scan_results/ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_2mt.txt")
  
show_gamma = False
use_log = False

@click.command()
@click.argument("plot_type")
def make_plot(plot_type):

  style()
    
  if plot_type == "scatter":
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    diff = fun_diff(data_mu_2, data_mu_05, data_mu_1, show_gamma=False, norm =False, sort=True, use_abs = False)
    vmax = 20
    size = 20
    title =  r'$T_C(Q=\frac{1}{2}m_t)-T_C(Q=2m_t)$'
    labels=[r'$M_s$ (GeV)', r'$\lambda_{S}$', r'$\lambda_{hs}$']
      
      
    def make_plot(ax, nx, ny, cbar=False):
      map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,4], vmax = 20, cmap="rainbow", edgecolor='none', s=size, alpha=1, rasterized=True)
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
        diff = fun_diff(data_mu_2, data_mu_05, data_mu_1, show_gamma=False, norm =False, use_abs = False)
        ax.scatter(diff[:,5], diff[:,4], c='seagreen', alpha=1, edgecolor='none', rasterized=True)
        ax.set_xlabel(r'$T_C(Q=m_t)$')
        ax.set_ylabel(r'$T_C(Q=\frac{1}{2}m_t)-T_C(Q=2m_t)$')
      elif nn == 1 :
        diff = fun_diff(data_mu_2, data_mu_05, data_mu_1, show_gamma=False, norm =True, use_abs = False)
        ax.scatter(diff[:,5], diff[:,4], c='seagreen', alpha=1, edgecolor='none', rasterized=True)
        ax.set_xlabel(r'$T_C(Q=m_t)$')
        ax.set_ylabel(r'$[T_C(Q=\frac{1}{2}m_t)-T_C(Q=2m_t)]/T_C(Q=m_t)$')
      else:
        diff = fun_diff(data_mu_2, data_mu_05, data_mu_1, show_gamma=True, norm =False, use_abs = False)
        ax.scatter(diff[:,5], diff[:,4], c='seagreen', alpha=1, edgecolor='none', rasterized=True)
        ax.set_xlabel(r'$\gamma(Q=m_t)$')
        ax.set_ylabel(r'$\gamma(Q=\frac{1}{2}m_t)-\gamma(Q=2m_t)$')
        ax.set_ylim(-5,0.1)
        ax.set_xlim(-0.1,7)
            
    make_plot(axs[0,0], 1, 2)
    make_plot(axs[0,1], 0, 1)
    make_plot(axs[0,2], 0, 2, cbar=True)
    no_color(axs[1,0], 0)
    no_color(axs[1,1], 1)
    no_color(axs[1,2], 2)
    plt.savefig('3d_mu_scatter.pdf')
  
  else:
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    
    show_gamma = plot_type == "gamma"
    def hist_plot(nn):
      show_gamma = (nn==1)
      ax = axs[nn,0]
      diff = fun_diff(data_mu_2, data_mu_05, data_mu_1, show_gamma=show_gamma, norm =False, use_abs = False)
      ax.grid(axis='y', alpha=0.75)
      ax.hist(x=diff[:,4], bins=30, color='seagreen', log=True, alpha=1, rwidth=0.85)
      if show_gamma:
        ax.set_xlabel(r'$\gamma(Q=\frac{1}{2}m_t)-\gamma(Q=2m_t)$')
      else:
        ax.set_xlabel(r'$T_C(Q=\frac{1}{2}m_t)-T_C(Q=2m_t)$')
      
      ax.set_ylabel("Number of samples")
    
      ax = axs[nn,1]
      diff = fun_diff(data_mu_2, data_mu_05, data_mu_1, show_gamma=show_gamma, norm =True, use_abs = False)
      ax.grid(axis='y', alpha=0.75)
      ax.hist(x=diff[:,4], bins=30, color='seagreen', log=True, alpha=1, rwidth=0.85)
      if show_gamma:
        ax.set_xlabel(r'$[\gamma(Q=\frac{1}{2}m_t)-\gamma(Q=2m_t)]/\gamma(Q=m_t)$')
      else:
        ax.set_xlabel(r'$[T_C(Q=\frac{1}{2}m_t)-T_C(Q=2m_t)]/T_C(Q=m_t)$')
      ax.set_ylabel("Number of samples")
    
    hist_plot(0)
    hist_plot(1)
    fig.tight_layout()
    plt.savefig('3d_mu_hist.pdf')
  
  plt.show()
 
if __name__ == "__main__":
    make_plot()
    

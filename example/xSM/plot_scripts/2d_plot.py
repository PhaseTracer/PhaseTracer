"""
2d plots for paper
==================
"""

import matplotlib.pyplot as plt
import numpy as np
from plot_fun import fun_diff
from style import style


style()
plot_xi = True
plot_scale = False
plot_scheme = False
plot_daisy = False


show_deltaT = False
show_deltagamma = False


if plot_xi:
  figure_name = "xi"
if plot_scale:
  figure_name = "scale"
if plot_scheme:
  figure_name = "scheme"

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
  
  show_data = fun_diff(data2, data1, data1, show_gamma=(plot_daisy and show_deltagamma), norm=(plot_daisy and show_deltagamma) )
    
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

figname = '2d_scan'
if plot_xi:
  figname += '_xi'
elif plot_scale:
  figname += '_scale'
elif plot_scheme:
  figname += '_scheme'
elif plot_daisy:
  figname += '_daisy'

if show_deltaT:
  if plot_daisy and show_deltagamma:
    figname += '_deltagamma'
  else:
    figname += '_deltaT'
else:
  figname += '_T'
  
plt.savefig(figname+'.pdf')

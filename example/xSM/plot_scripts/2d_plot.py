import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from plot_fun import fun_gamma, fun_diff, loaddata
  
plot_xi = True
plot_scale = not plot_xi

show_deltaT = True

if plot_xi:
  figure_name = "xi"
if plot_scale:
  figure_name = "scale"



if show_deltaT:
  if plot_xi:
    title=r"$|T_C^{(\xi=3)}-T_C^{(\xi=0)}|$ [GeV]"
  if plot_scale:
    title=r"$|T_C^{(Q=\frac{1}{2}m_t)}-T_C^{(Q=2m_t)}|$ [GeV]"
  cm = 'rainbow'
else:
  title=r"$T_C^{(\xi=0)}$ [GeV]"
  cm = 'cool'

labels=[r'$M_s$ [GeV]', r'$\lambda_{S}$', r'$\lambda_{HS}$']

def make_plot(ax, par):
  #############################
  if plot_xi:
    data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
    data2 = np.loadtxt("../2d_scan/"+par+"_xi3.txt")
    vmax = 10
  if plot_scale:
    data1 = np.loadtxt("../2d_scan/"+par+"_05mt.txt")
    data2 = np.loadtxt("../2d_scan/"+par+"_2mt.txt")
    vmax = 15
    
  show_data = fun_diff(data2, data1, data1)
    
  # 0=lambda_hs, 1=lambda_s, 2=ms
  if par=="lhs_ls":    
    nx = 1
    ny = 2
    label = r"$M_s$=62.6 GeV"
    
  if par=="ms_lhs":
    nx = 0
    ny = 2
    label = r"$\lambda_{S}$=0.1"

  if par=="ms_ls":
    nx = 0
    ny = 1
    label = r"$\lambda_{HS}=0.3$"

  xmin = min(show_data[:,nx])
  xmax = max(show_data[:,nx])
  ymin = min(show_data[:,ny])
  ymax = max(show_data[:,ny])
    
  if show_deltaT:
    map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=abs(show_data[:,4]), cmap=cm, edgecolor='none', s=5, vmax=vmax, alpha=1)
  else:
    map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=show_data[:,5], cmap=cm, s=2, vmax=150, alpha=1)
    
  ax.set_xlabel(labels[nx])
  ax.set_ylabel(labels[ny])
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)
  ax.set_title(title)
  ax.text(xmin+0.1*(xmax-xmin),ymax-0.1*(ymax-ymin),label)

  clb = plt.colorbar(map1, ax=ax)


fig, axs = plt.subplots(1,3, figsize=(15, 5))

make_plot(axs[0], par="lhs_ls")
make_plot(axs[1], par="ms_ls")
make_plot(axs[2], par="ms_lhs")

fig.tight_layout()

#plt.subplots_adjust(hspace=0.2)

figname = '2d_scan'
if plot_xi:
  figname += '_xi'
elif plot_scale:
  figname += '_scale'
if show_deltaT:
  figname += '_deltaT'
else:
  figname += '_T'
  
plt.savefig(figname+'.pdf')
plt.show()

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from plot_fun import fun_gamma, fun_diff, loaddata
  
show_gamma = False

  
figure_name = "xi"
show_deltaT = False




if show_deltaT:
  title=r"$|T_C^{(\xi=3)}-T_C^{(\xi=0)}|$ [GeV]"
  cm = 'rainbow'
else:
  title=r"$T_C^{(\xi=0)}$ [GeV]"
  cm = 'cool'

labels=[r'$M_s$ [GeV]', r'$\lambda_{S}$', r'$\lambda_{HS}$']

def make_plot(ax, par):
  #############################
  data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
  data2 = np.loadtxt("../2d_scan/"+par+"_xi3.txt")
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
    map1 = ax.scatter(show_data[:,nx], show_data[:,ny], c=abs(show_data[:,4]), cmap=cm, s=2, vmax=10, alpha=1)
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

figname = '2d_scan_xi'
if show_deltaT:
  figname += '_deltaT'
else:
  figname += '_T'
  
plt.savefig(figname+'.pdf')
plt.show()

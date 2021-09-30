import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
sys.path.append( '../' )
from plot_fun import fun_gamma, fun_diff, loaddata

cmap = cm.get_cmap('rainbow')

fig, axs = plt.subplots(2, 3, figsize=(10, 6))

def plot_for_1d(data, x_num, label, column):

  sel = data[:,3]>0

  x = data[:,x_num]
  TC = data[:,4]
  gamma = fun_gamma(data)

  ax = axs[0,column]
  ax.plot(x[sel], TC[sel], label=label, alpha=1)

  ax = axs[1,column]
  ax.plot(x[sel], gamma[sel], label=label, alpha=1)

names = [ ["default", r"MS, $Q=m_t$"],
          ["2mt", r"MS, $Q=2m_t$"],
          ["05mt", r"MS, $Q=m_t/2$"],
          ["PRM", "PRM $Q=m_t$"],
          ["PRM_2mt", "PRM $Q=2m_t$"],
          ["PRM_05mt", "PRM $Q=m_t/2$"]]

for name in names:
  plot_for_1d(np.loadtxt("../1d_bks/lambda_hs_"+name[0]+".txt"), 2, name[1], 0)
  plot_for_1d(np.loadtxt("../1d_bks/lambda_s_"+name[0]+".txt"), 1, name[1], 1)
  plot_for_1d(np.loadtxt("../1d_bks/m_s_"+name[0]+".txt"), 0, name[1], 2)
  

for ii in range(2):
  for jj in range(3):
    axs[ii,jj].grid(axis='x', alpha=0.75)
    axs[ii,jj].grid(axis='y', alpha=0.75)
    
    if ii == 0:
      axs[ii,jj].set_ylabel(r"$T_C$ (GeV)")
      if jj<1:
        axs[ii,jj].legend(loc=3)
      else:
        axs[ii,jj].legend(loc=4)
    else:
      axs[ii,jj].set_ylim(0,10)
      axs[ii,jj].set_ylabel(r"$\gamma_{\rm EW}$")
      if jj<1:
        axs[ii,jj].legend(loc=2)
      else:
        axs[ii,jj].legend(loc=1)
        
    if jj == 0:
      axs[ii,jj].set_xlabel(r"$\lambda_{hs}$")
    elif jj == 1:
      axs[ii,jj].set_xlabel(r"$\lambda_{s}$")
    else:
      axs[ii,jj].set_xlabel(r"$m_{s}$ (GeV)")

fig.tight_layout()
plt.savefig('1d_bks.png')
plt.show()

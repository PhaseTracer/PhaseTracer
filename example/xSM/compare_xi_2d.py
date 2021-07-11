import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_fun import fun_gamma, fun_diff, loaddata

if False:
  scheme = "OSlike_2d"
  if sys.argv[1] == '1':
    data_xi_0 = np.loadtxt("onshell/xSM_OSlike_xi_0_ls_fixed.txt")
    data_xi_1 = np.loadtxt("onshell/xSM_OSlike_xi_1_ls_fixed.txt")
    data_xi_3 = np.loadtxt("onshell/xSM_OSlike_xi_3_ls_fixed.txt")
  elif sys.argv[1] == '3':
    data_xi_0 = np.loadtxt("onshell/xSM_OSlike_xi_0_lhs_fixed.txt")
    data_xi_1 = np.loadtxt("onshell/xSM_OSlike_xi_1_lhs_fixed.txt")
    data_xi_3 = np.loadtxt("onshell/xSM_OSlike_xi_3_lhs_fixed.txt")
  else:
    data_xi_0 = np.loadtxt("onshell/xSM_OSlike_xi_0_ms_fixed.txt")
    data_xi_1 = np.loadtxt("onshell/xSM_OSlike_xi_1_ms_fixed.txt")
    data_xi_3 = np.loadtxt("onshell/xSM_OSlike_xi_3_ms_fixed.txt")

  marker_size = 100
else:
  scheme = "MSbar"
  if sys.argv[1] == '1':
    data_xi_0 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_0_ls_fixed.txt")
    data_xi_1 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_1_ls_fixed.txt")
    data_xi_3 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_3_ls_fixed.txt")
  elif sys.argv[1] == '3':
    data_xi_0 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_0_lhs_fixed.txt")
    data_xi_1 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_1_lhs_fixed.txt")
    data_xi_3 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_3_lhs_fixed.txt")
  else:
    data_xi_0 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_0_ms_fixed.txt")
    data_xi_1 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_1_ms_fixed.txt")
    data_xi_3 = np.loadtxt("gauge_dependence/xSM_MSbar_xi_3_ms_fixed.txt")
  marker_size = 1





show_gamma = False
use_log = False

if len(sys.argv)<2 :
  par = "lhs_ls"
else:
  if sys.argv[1] == '1':
    par = "lhs_ms"
  elif sys.argv[1] == '2':
    par = "lhs_ls"
  elif sys.argv[1] == '3':
    par = "ls_ms"

if par=="lhs_ls":    
    nx = 1
    ny = 2
    nz = 0
    xlabel = r'$\lambda_{S}$'
    ylabel = r'$\lambda_{HS}$'
    zlabel = r'$m_{S}=65$ GeV'
    xmin = 0
    xmax = 0.3
    ymin = 0.1
    ymax = 0.5
    figure_name = "lhs_ls"

if par=="lhs_ms":    
    nx = 0
    ny = 2
    nz = 1
    xlabel = r'$m_{S}$'
    ylabel = r'$\lambda_{HS}$'
    zlabel = r'$\lambda_{S}=0.1$'
    xmin = 10
    xmax = 110
    ymin = 0.1
    ymax = 0.5
    figure_name = "lhs_ms"

if par=="ls_ms":    
    nx = 0
    ny = 1
    nz = 2
    xlabel = r'$m_{S}$'
    ylabel = r'$\lambda_{S}$'
    zlabel = r'$\lambda_{HS}=0.3$'
    xmin = 10
    xmax = 110
    ymin = 0
    ymax = 0.3
    figure_name = "ls_ms"

def get_griddata(px,py,nx,ny,c2):

    dx = (max(px) - min(px))/nx/2.0
    dy = (max(py) - min(py))/ny/2.0

    xi = np.linspace(min(px-dx), max(px+dx), nx)
    yi = np.linspace(min(py-dy), max(py+dy), ny)

    xf = np.zeros((nx,ny))
    yf = np.zeros((nx,ny))
    zf = np.zeros((nx,ny))

    for ii,ix in enumerate(xi):
        for jj,jy in enumerate(yi):
            xf[ii,jj] = ix
            yf[ii,jj] = jy
            if len(c2[(px<ix+dx)&(px>ix-dx)&(py<jy+dy)&(py>jy-dy)])>0:
                value = max(c2[(px<ix+dx)&(px>ix-dx)&(py<jy+dy)&(py>jy-dy)])
            else:
                value = -10
            zf[ii,jj] = value

    return [xf,yf,zf]


if True:

  diff = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma, norm =False, use_abs=False)
  if use_log:
    diff[:,3] = np.log10(diff[:,3]+1E-10) 
    vmin = -5
    vmax = -1
  else:
    vmin = 0
    vmax = 0.2
  if use_log:
    diff[:,3] = np.log10(diff[:,3]+1E-10) 
    vmin = -5
    vmax = -1
  else:
    vmin = 0
    vmax = 0.2

  fig, axs = plt.subplots(2, 1, figsize=(5, 9))

  ax = axs[0]
  map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,3], cmap="rainbow", s=marker_size, marker=".", alpha=1)

  ax.text(xmin + 0.1*(xmax-xmin), ymax - 0.1*(ymax-ymin), zlabel, fontsize=12, color='k')

  clb = plt.colorbar(map, ax=ax)


  label = (r'$T_C^{(\xi=0)}-T_C^{(\xi=3)}$' )
  if use_log:
    label = r'log$_{10}$('+label+')'

  ax.set_title(label=label)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)

  ax = axs[1]
#  map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,4], cmap="cool", s=marker_size, marker=".", alpha=1)
  sel = (abs(data_xi_0[:,5])>1) & (abs(data_xi_0[:,6])<1) & (abs(data_xi_0[:,7])<1) & (abs(data_xi_0[:,8])>1)
  map = ax.scatter(data_xi_0[sel][:,nx], data_xi_0[sel][:,ny], c=data_xi_0[sel][:,4], cmap="cool", s=marker_size, marker=".", alpha=1)
  
  clb = plt.colorbar(map, ax=ax)

  ax.set_title(label=(r"$T_C^{(\xi=0)}$") )
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)

  ax.text(xmin + 0.1*(xmax-xmin), ymax - 0.1*(ymax-ymin), zlabel, fontsize=12, color='k')
  
  fig.tight_layout()
  plt.savefig(scheme+'_gauge_'+figure_name + ('_gamma' if show_gamma else '_Tc') + ( '_log' if use_log else '') + '.png')
  

plt.show()

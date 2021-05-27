import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_fun import fun_gamma, fun_diff, loaddata


show_gamma = False
use_log = False

#data_AE = np.loadtxt("scan_results/xSM_MSbar_daisy_ArnoldEspinosa.txt")
#data_PW = np.loadtxt("scan_results/xSM_MSbar_daisy_Parwani.txt")
#diff = fun_diff(data_AE, data_PW, data_PW, show_gamma)

data_xi_0 = np.loadtxt("random_scan_results/xSM_MSbar_xi_0.txt")
data_xi_01 = np.loadtxt("random_scan_results/xSM_MSbar_xi_01.txt")
data_xi_1 = np.loadtxt("random_scan_results/xSM_MSbar_xi_1.txt")
data_xi_3 = np.loadtxt("random_scan_results/xSM_MSbar_xi_3.txt")
diff = fun_diff(data_xi_01, data_xi_0, data_xi_1, show_gamma)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))


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
    xlabel = r'$\lambda_{S}$'
    ylabel = r'$\lambda_{HS}$'
    zlabel = r'$m_{S}$'
    xmin = 0
    xmax = 0.3
    ymin = 0.1
    ymax = 0.5
    figure_name = "lhs_ls"

if par=="lhs_ms":    
    nx = 0
    ny = 2
    xlabel = r'$m_{S}$'
    ylabel = r'$\lambda_{HS}$'
    zlabel = r'$\lambda_{S}$'
    xmin = 10
    xmax = 110
    ymin = 0.1
    ymax = 0.5
    figure_name = "lhs_ms"

if par=="ls_ms":    
    nx = 0
    ny = 1
    xlabel = r'$m_{S}$'
    ylabel = r'$\lambda_{S}$'
    zlabel = r'$\lambda_{HS}$'
    xmin = 10
    xmax = 110
    ymin = 0
    ymax = 0.3
    figure_name = "ls_ms"


if use_log:
  diff[:,3] = np.log10(diff[:,3]+1E-10) 
  vmin = -5
  vmax = -1
else:
  vmin = 0
  vmax = 0.2


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

ax = axs[0]
map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,3], cmap="autumn", s=20, marker="s", vmin=vmin, vmax =0.1,alpha=1)

#[xf,yf,zf] = get_griddata(diff[:,nx], diff[:,ny], 50, 50, diff[:,3])
#map = ax.scatter(xf[zf>-9],yf[zf>-9], c=zf[zf>-9], cmap="autumn", s=20, marker="s", vmin=vmin, vmax =0.1, alpha=1)

clb = plt.colorbar(map, ax=ax)

label = (r'$|\Delta \gamma_{\rm EW}| / \gamma_{\rm EW}^{bk}$' if show_gamma else r'$|\Delta T_C|/T_C$' )
if use_log:
  label = r'log$_{10}$('+label+')'

ax.set_title(label=label)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)


ax = axs[1]
map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,0], cmap="summer", s=20, marker="s", alpha=1)
clb = plt.colorbar(map, ax=ax)

ax.set_title(label="Corresponding "+zlabel)
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax = axs[2]
map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,4], cmap="winter", s=20, marker="s", alpha=1)
clb = plt.colorbar(map, ax=ax)

ax.set_title(label="Corresponding "+ (r"$\gamma_{\rm EW}^{bk}$" if show_gamma else "$T_C^{BK}$") )
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)


#fig.tight_layout()

plt.savefig('compare_'+figure_name + ('_gamma' if show_gamma else '_Tc') + ( '_log' if use_log else '') + '.png')
plt.show()

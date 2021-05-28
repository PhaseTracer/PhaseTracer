import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_fun import fun_gamma, fun_diff, loaddata


fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))

for ax in axs:
  ax.grid(axis='x', alpha=0.75)
  ax.grid(axis='y', alpha=0.75)

def plot_for_gauge(data_MS, color):
  xi = data_MS[:,10]
  TC_MS = data_MS[:,4]
  gamma_MS = fun_gamma(data_MS)

  ax = axs[0]
  ax.plot(xi, TC_MS,  c=color, alpha=1)

  ax = axs[1]
  ax.plot(xi, gamma_MS,  c=color, alpha=1)


plot_for_gauge(np.loadtxt("gauge_dependence/MSbar_BK1.txt"),  "orange")
plot_for_gauge(np.loadtxt("gauge_dependence/MSbar_BK2.txt"),  "r")
plot_for_gauge(np.loadtxt("gauge_dependence/MSbar_BK3.txt"),  "purple")
plot_for_gauge(np.loadtxt("gauge_dependence/MSbar_BK4.txt"),  "g")
plot_for_gauge(np.loadtxt("gauge_dependence/MSbar_BK5.txt"),  "b")
#plot_for_gauge(np.loadtxt("gauge_dependence/MSbar_62.5_0.11_0.25.txt"))


axs[0].set_xlabel(r"$\xi$")
axs[0].set_ylabel(r"$T_C$")

axs[1].set_xlabel(r"$\xi$")
axs[1].set_ylabel(r"$\gamma$")

fig.tight_layout()
plt.savefig('gauge_1d.png')

#label = (r'$|\Delta \gamma_{\rm EW}| / \gamma_{\rm EW}^{bk}$' if show_gamma else r'$|\Delta T_C|/T_C$' )
#if use_log:
#  label = r'log$_{10}$('+label+')'

#ax.set_title(label=label)
#ax.set_xlabel(xlabel)
#ax.set_ylabel(ylabel)
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)


#if len(sys.argv)<2 :
#  par = "lhs_ls"
#else:
#  if sys.argv[1] == '1':
#    par = "lhs_ms"
#  elif sys.argv[1] == '2':
#    par = "lhs_ls"
#  elif sys.argv[1] == '3':
#    par = "ls_ms"

#if par=="lhs_ls":    
#    nx = 1
#    ny = 2
#    xlabel = r'$\lambda_{S}$'
#    ylabel = r'$\lambda_{HS}$'
#    zlabel = r'$m_{S}$'
#    xmin = 0
#    xmax = 0.3
#    ymin = 0.1
#    ymax = 0.5
#    figure_name = "lhs_ls"

#if par=="lhs_ms":    
#    nx = 0
#    ny = 2
#    xlabel = r'$m_{S}$'
#    ylabel = r'$\lambda_{HS}$'
#    zlabel = r'$\lambda_{S}$'
#    xmin = 10
#    xmax = 110
#    ymin = 0.1
#    ymax = 0.5
#    figure_name = "lhs_ms"

#if par=="ls_ms":    
#    nx = 0
#    ny = 1
#    xlabel = r'$m_{S}$'
#    ylabel = r'$\lambda_{S}$'
#    zlabel = r'$\lambda_{HS}$'
#    xmin = 10
#    xmax = 110
#    ymin = 0
#    ymax = 0.3
#    figure_name = "ls_ms"


#if use_log:
#  diff[:,3] = np.log10(diff[:,3]+1E-10) 
#  vmin = -5
#  vmax = -1
#else:
#  vmin = 0
#  vmax = 0.2


#def get_griddata(px,py,nx,ny,c2):

#    dx = (max(px) - min(px))/nx/2.0
#    dy = (max(py) - min(py))/ny/2.0

#    xi = np.linspace(min(px-dx), max(px+dx), nx)
#    yi = np.linspace(min(py-dy), max(py+dy), ny)

#    xf = np.zeros((nx,ny))
#    yf = np.zeros((nx,ny))
#    zf = np.zeros((nx,ny))

#    for ii,ix in enumerate(xi):
#        for jj,jy in enumerate(yi):
#            xf[ii,jj] = ix
#            yf[ii,jj] = jy
#            if len(c2[(px<ix+dx)&(px>ix-dx)&(py<jy+dy)&(py>jy-dy)])>0:
#                value = max(c2[(px<ix+dx)&(px>ix-dx)&(py<jy+dy)&(py>jy-dy)])
#            else:
#                value = -10
#            zf[ii,jj] = value

#    return [xf,yf,zf]

#ax = axs[0]
#map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,3], cmap="autumn", s=20, marker="s", vmin=vmin, vmax =0.1,alpha=1)

##[xf,yf,zf] = get_griddata(diff[:,nx], diff[:,ny], 50, 50, diff[:,3])
##map = ax.scatter(xf[zf>-9],yf[zf>-9], c=zf[zf>-9], cmap="autumn", s=20, marker="s", vmin=vmin, vmax =0.1, alpha=1)

#clb = plt.colorbar(map, ax=ax)

#label = (r'$|\Delta \gamma_{\rm EW}| / \gamma_{\rm EW}^{bk}$' if show_gamma else r'$|\Delta T_C|/T_C$' )
#if use_log:
#  label = r'log$_{10}$('+label+')'

#ax.set_title(label=label)
#ax.set_xlabel(xlabel)
#ax.set_ylabel(ylabel)
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)


#ax = axs[1]
#map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,0], cmap="summer", s=20, marker="s", alpha=1)
#clb = plt.colorbar(map, ax=ax)

#ax.set_title(label="Corresponding "+zlabel)
#ax.set_xlabel(xlabel)
#ax.set_ylabel(ylabel)
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)

#ax = axs[2]
#map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,4], cmap="winter", s=20, marker="s", alpha=1)
#clb = plt.colorbar(map, ax=ax)

#ax.set_title(label="Corresponding "+ (r"$\gamma_{\rm EW}^{bk}$" if show_gamma else "$T_C^{BK}$") )
#ax.set_xlabel(xlabel)
#ax.set_ylabel(ylabel)
#ax.set_xlim(xmin,xmax)
#ax.set_ylim(ymin,ymax)


#fig.tight_layout()

#plt.savefig('compare_'+figure_name + ('_gamma' if show_gamma else '_Tc') + ( '_log' if use_log else '') + '.png')
plt.show()

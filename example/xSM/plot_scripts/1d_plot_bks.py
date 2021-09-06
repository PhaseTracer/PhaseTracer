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

names = [ ["default", "default"],
          ["2mt", r"$Q=2m_t$"],
          ["OSlike", r"OS-like"],
          ["xi3", r"$\xi=3$"],
          ["HT", r"HT"],
          ["PRM","PRM"]]

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
      axs[ii,jj].set_xlabel(r"$m_{s}$")

fig.tight_layout()
plt.savefig('1d_bks.png')

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

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_fun import fun_gamma, fun_diff, loaddata


data_xi_0 = np.loadtxt("random_scan_results/xSM_MSbar_xi_0.txt")
data_xi_1 = np.loadtxt("random_scan_results/xSM_MSbar_xi_1.txt")
data_xi_3 = np.loadtxt("random_scan_results/xSM_MSbar_xi_3.txt")
  
#data_xi_1 = np.loadtxt("random_scan_results_new/xSM_MSbar_xi_1_use_1L_EWSB_in_0L_mass.txt")
#data_xi_3 = np.loadtxt("random_scan_results_new/xSM_MSbar_xi_3_use_1L_EWSB_in_0L_mass.txt")



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
    zlabel = r'$m_{S}$'
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
    zlabel = r'$\lambda_{S}$'
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
    zlabel = r'$\lambda_{HS}$'
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

  diff = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma)
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

  fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))

  ax = axs[0]
  map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,3], cmap="rainbow", vmax = 0.1, s=30, marker=".", vmin=vmin,alpha=1)

  #[xf,yf,zf] = get_griddata(diff[:,nx], diff[:,ny], 50, 50, diff[:,3])
  #map = ax.scatter(xf[zf>-9],yf[zf>-9], c=zf[zf>-9], cmap="autumn", s=20, marker="s", vmin=vmin, vmax =0.1, alpha=1)

  clb = plt.colorbar(map, ax=ax)


  label = (r'$|\gamma^{(\xi=0)} - \gamma^{(\xi=3)}| / \gamma^{(\xi=1)}$' if show_gamma else r'$|T_C^{(\xi=0)}-T_C^{(\xi=3)}|/T_C^{(\xi=1)}$' )
  if use_log:
    label = r'log$_{10}$('+label+')'

  ax.set_title(label=label)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)


  ax = axs[1]
  map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,nz], cmap="autumn", s=30, marker=".", alpha=1)
  clb = plt.colorbar(map, ax=ax)

  ax.set_title(label="Corresponding "+zlabel)
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)

  ax = axs[2]
  map = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,4], cmap="winter", s=30, marker=".", alpha=1)
  clb = plt.colorbar(map, ax=ax)

  ax.set_title(label="Corresponding "+ (r"$\gamma^{(\xi=1)}$" if show_gamma else r"$T_C^{(\xi=1)}$") )
  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)

  fig.tight_layout()
  plt.savefig('gauge_'+figure_name + ('_gamma' if show_gamma else '_Tc') + ( '_log' if use_log else '') + '.png')


if False:

  norm = True
  ZoomIn = False

  label_TC = r'$[T_C^{(\xi=0)}-T_C^{(\xi=3)}]/T_C^{(\xi=1)}$' if norm else r'$T_C^{(\xi=0)} - T_C^{(\xi=3)}$'
  label_gamma = r'$[\gamma^{(\xi=0)}-\gamma^{(\xi=3)}]/\gamma^{(\xi=1)}$' if norm else r'$\gamma^{(\xi=0)} - \gamma^{(\xi=3)}$'
  
  diff_for_gamma = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma=True, norm = norm, use_abs = False, sort=False)
  diff_gamma = diff_for_gamma[:,3]
  diff_for_TC = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma=False, norm = norm, use_abs = False, sort=False)
  diff_TC = diff_for_TC[:,3]
  
  fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))

  ax = axs[0]
  ax.grid(axis='y', alpha=0.75)
  
  ax.hist(x=diff_TC, bins=50, color='g', log=True, alpha=0.7, rwidth=0.85)

  ax.set_xlabel(label_TC)
  ax.set_ylabel("Number of samples")

  ax = axs[1]
  ax.grid(axis='y', alpha=0.75)
  if ZoomIn:
    ax.hist(x=diff_gamma, bins=50, range = [-1,1], color='g', log=True, alpha=0.7, rwidth=0.85)
  else:
    ax.hist(x=diff_gamma, bins=50, color='g', log=True, alpha=0.7, rwidth=0.85)
  ax.set_xlabel(label_gamma)
  ax.set_ylabel("Number of samples")
  
  
  ax = axs[2]
  ax.grid(axis='x', alpha=0.75)
  ax.grid(axis='y', alpha=0.75)
  ax.scatter(diff_TC, diff_gamma, c='g', alpha=0.7)
  ax.set_xlabel(label_TC)
  ax.set_ylabel(label_gamma)
  if ZoomIn:
    ax.set_xlim(-0.1,0.15)
    ax.set_ylim(-1,2)

  fig.tight_layout()
  plt.savefig('gaugue_hist'+ ('_norm' if norm else '') + ('_ZoomIn' if ZoomIn else '') + '.png')
  
  
if False:

  norm = True
  ZoomIn = False

  one_point = False

  label_TC = r'$[T_C^{(\xi=0)}-T_C^{(\xi=3)}]/T_C^{(\xi=1)}$' if norm else r'$T_C^{(\xi=0)} - T_C^{(\xi=3)}$'
  label_gamma = r'$[\gamma^{(\xi=0)}-\gamma^{(\xi=3)}]/\gamma^{(\xi=1)}$' if norm else r'$\gamma^{(\xi=0)} - \gamma^{(\xi=3)}$'
  


  diff_for_gamma = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma=True, norm = norm, use_abs = False, sort=False)
  diff_gamma = diff_for_gamma[:,3]
  diff_for_TC = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma=False, norm = norm, use_abs = False, sort=False)
  diff_TC = diff_for_TC[:,3]

  fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))

  if one_point:
    sel_1 = diff_TC < - 0.4
    sel_2 = abs(diff_TC+0.2) < 1E-3
    sel_3 = (abs(diff_TC+0.1) < 1E-3) & (abs(diff_gamma-0.1) < 1E-2)
    sel_4 = (abs(diff_TC+0.06) < 5E-4) & (abs(diff_gamma-0.1) < 5E-3)
    sel_5 = (diff_TC > -0.048) & (abs(diff_gamma-0.15) < 1.2E-3)
    sel_6 = diff_TC > -0.0135
  else:
    sel_1 = diff_TC < - 0.21
    sel_2 = (diff_TC > - 0.21) & (diff_TC < - 0.15)
    sel_3 = (diff_TC > - 0.15) & (diff_TC < - 0.05) & (diff_gamma/diff_TC > -1.1)
    sel_4 =  (diff_gamma/diff_TC < -1.3) & (diff_gamma < -6*diff_TC - 0.3) 
    sel_5 = (diff_TC < - 0.04) & (diff_gamma > -6*diff_TC - 0.15) 
    sel_6 = diff_TC > -0.02




  print "BK = 'BK1'"
  print "ms = ",        diff_for_gamma[:,0][sel_1]
  print "lambda_s = ",  diff_for_gamma[:,1][sel_1]
  print "lambda_hs = ", diff_for_gamma[:,2][sel_1]
  if len(diff_for_gamma[:,0][sel_1])>1: print diff_for_gamma[sel_1]
  
  print "BK = 'BK2'"
  print "ms = ",        diff_for_gamma[:,0][sel_2]
  print "lambda_s = ",  diff_for_gamma[:,1][sel_2]
  print "lambda_hs = ", diff_for_gamma[:,2][sel_2]
  if len(diff_for_gamma[:,0][sel_2])>1: print diff_for_gamma[sel_2]
  
  print "BK = 'BK3'"
  print "ms = ",        diff_for_gamma[:,0][sel_3]
  print "lambda_s = ",  diff_for_gamma[:,1][sel_3]
  print "lambda_hs = ", diff_for_gamma[:,2][sel_3]
  if len(diff_for_gamma[:,0][sel_3])>1: print diff_for_gamma[sel_3]
  
  print "BK = 'BK4'"
  print "ms = ",        diff_for_gamma[:,0][sel_4]
  print "lambda_s = ",  diff_for_gamma[:,1][sel_4]
  print "lambda_hs = ", diff_for_gamma[:,2][sel_4]
  if len(diff_for_gamma[:,0][sel_4])>1: print diff_for_gamma[sel_4]
    
  print "BK = 'BK5'"
  print "ms = ",        diff_for_gamma[:,0][sel_5]
  print "lambda_s = ",  diff_for_gamma[:,1][sel_5]
  print "lambda_hs = ", diff_for_gamma[:,2][sel_5]
  if len(diff_for_gamma[:,0][sel_5])>1: print diff_for_gamma[sel_5]

  print "BK = 'BK6'"
  print "ms = ",        diff_for_gamma[:,0][sel_6]
  print "lambda_s = ",  diff_for_gamma[:,1][sel_6]
  print "lambda_hs = ", diff_for_gamma[:,2][sel_6]
  if len(diff_for_TC[:,0][sel_6])>1: print diff_for_TC[sel_6]


  norm = False
  diff_for_gamma = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma=True, norm = norm, use_abs = False, sort=False)
  diff_gamma = diff_for_gamma[:,3]
  diff_for_TC = fun_diff(data_xi_0, data_xi_3, data_xi_1, show_gamma=False, norm = norm, use_abs = False, sort=False)
  diff_TC = diff_for_TC[:,3]
  label_TC = r'$[T_C^{(\xi=0)}-T_C^{(\xi=3)}]/T_C^{(\xi=1)}$' if norm else r'$T_C^{(\xi=0)} - T_C^{(\xi=3)}$'
  label_gamma = r'$[\gamma^{(\xi=0)}-\gamma^{(\xi=3)}]/\gamma^{(\xi=1)}$' if norm else r'$\gamma^{(\xi=0)} - \gamma^{(\xi=3)}$'
  
#  bk6 = (data_xi_0[:,0] == diff_for_gamma[:,0][sel_6] ) & (data_xi_0[:,1] == diff_for_gamma[:,1][sel_6] )
#  
#  print data_xi_0[bk6]
#  print data_xi_1[bk6]
#  print data_xi_3[bk6]

#  print diff_for_TC[sel_6]  
#  print diff_for_gamma[sel_6]

#  
#  print diff_TC[sel_6], diff_gamma[sel_6]
#  

  ax = axs[0]
  ax.grid(axis='x', alpha=0.75)
  ax.grid(axis='y', alpha=0.75)
  ax.scatter(diff_for_gamma[:,1], diff_for_gamma[:,2], c='gray', alpha=0.2)
  
  ax.scatter(diff_for_gamma[:,1][sel_3], diff_for_gamma[:,2][sel_3], c='g', alpha=1)
  ax.scatter(diff_for_gamma[:,1][sel_5], diff_for_gamma[:,2][sel_5], c='b', alpha=1)
  ax.scatter(diff_for_gamma[:,1][sel_6], diff_for_gamma[:,2][sel_6], c='purple', alpha=1)
  ax.scatter(diff_for_gamma[:,1][sel_4], diff_for_gamma[:,2][sel_4], c='cyan', alpha=1)
  
  ax.scatter(diff_for_gamma[:,1][sel_2], diff_for_gamma[:,2][sel_2], c='orange', alpha=1)
  
  ax.scatter(diff_for_gamma[:,1][sel_1], diff_for_gamma[:,2][sel_1], c='r', alpha=1)
  ax.set_xlabel(r'$\lambda_{S}$')
  ax.set_ylabel(r'$\lambda_{HS}$')


  ax = axs[1]
  ax.grid(axis='x', alpha=0.75)
  ax.grid(axis='y', alpha=0.75)
  ax.scatter(diff_for_gamma[:,0], diff_for_gamma[:,1], c='gray', alpha=0.2)

  
  ax.scatter(diff_for_gamma[:,0][sel_3], diff_for_gamma[:,1][sel_3], c='g', alpha=1)
  
  ax.scatter(diff_for_gamma[:,0][sel_5], diff_for_gamma[:,1][sel_5], c='b', alpha=1)
  ax.scatter(diff_for_gamma[:,0][sel_6], diff_for_gamma[:,1][sel_6], c='purple', alpha=1)
  
  ax.scatter(diff_for_gamma[:,0][sel_4], diff_for_gamma[:,1][sel_4], c='cyan', alpha=1)
  ax.scatter(diff_for_gamma[:,0][sel_2], diff_for_gamma[:,1][sel_2], c='orange', alpha=1)
  ax.scatter(diff_for_gamma[:,0][sel_1], diff_for_gamma[:,1][sel_1], c='r', alpha=1)
  
  ax.set_xlabel(r'$m_{S}$')
  ax.set_ylabel(r'$\lambda_{S}$')

    
  ax = axs[2]
  ax.grid(axis='x', alpha=0.75)
  ax.grid(axis='y', alpha=0.75)
  ax.scatter(diff_TC, diff_gamma, c='gray', alpha=0.2)
  
  ax.scatter(diff_TC[sel_1], diff_gamma[sel_1], c='r', alpha=1)
  ax.scatter(diff_TC[sel_2], diff_gamma[sel_2], c='orange', alpha=1)
  ax.scatter(diff_TC[sel_3], diff_gamma[sel_3], c='g', alpha=1)
  ax.scatter(diff_TC[sel_4], diff_gamma[sel_4], c='cyan', alpha=1)
  ax.scatter(diff_TC[sel_5], diff_gamma[sel_5], c='b', alpha=1)
  ax.scatter(diff_TC[sel_6], diff_gamma[sel_6], c='purple', alpha=1)
  
  
  ax.set_xlabel(label_TC)
  ax.set_ylabel(label_gamma)
  if ZoomIn:
    ax.set_ylim(-1,1)

  fig.tight_layout()
  plt.savefig('gaugue_sel_bk_points.png' if one_point else 'gaugue_classify.png')
  

plt.show()

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from useful_fun import fun_gamma, fun_diff, loaddata


if len(sys.argv)<3 :
  print "Use compare.py 0 1"
  sys.exit()

show_gamma = sys.argv[1] == '1'

cm = 'autumn_r' if show_gamma else 'autumn'

if sys.argv[2] == '1':
  compare_label = r'$Q=2m_t$'
  compare_name = "MSbar_ArnoldEspinosa_2mt_NoTreeEWSB"
  figure_name = "2mt"
elif sys.argv[2] == '2':
  compare_label = r'Goldstone resummation'
  compare_name = "MSbar_ArnoldEspinosa_mt_TreeEWSB"
  figure_name = "Goldstone"
elif sys.argv[2] == '3':
  compare_label = r'Parwani resummation'
  compare_name = "MSbar_Parwani_mt_NoTreeEWSB"
  figure_name = "Parwani"
elif sys.argv[2] == '4':
  compare_label = r'$\xi=1$'
  compare_name = "MSbar_ArnoldEspinosa_mt_NoTreeEWSB_xi"
  figure_name = "xi"

#label=r'On-shell like')
#label=r'High-$T$ approximation')
#label=r'PRM')
  
suptitle = "BK scenario $vs.$ ["+compare_label+"]"
print compare_label

data1 = np.loadtxt("Benchmarked/MSbar_ArnoldEspinosa_mt_NoTreeEWSB_lhs_ls.txt")
data2 = np.loadtxt("Compare/"+ compare_name +"_lhs_ls.txt")
diff, d1, d2 = fun_diff(data1, data2, show_gamma)
vmin = max(min(diff[:,3]),-0.1)
vmax = min(max(diff[:,3]),0.1)

def make_plot(ax, par):
  #############################
  data1 = np.loadtxt("Benchmarked/MSbar_ArnoldEspinosa_mt_NoTreeEWSB_"+par+".txt")
  data2 = np.loadtxt("Compare/"+ compare_name +"_"+par+".txt")

  diff, d1, d2 = fun_diff(data1, data2, show_gamma)
  
  # 0=lambda_hs, 1=lambda_s, 2=ms
  if par=="lhs_ls":    
    nx = 1
    ny = 0
    xlabel = r'$\lambda_{S}$'
    ylabel = r'$\lambda_{HS}$'
    xmin = 0
    xmax = 0.2
    ymin = 0.1
    ymax = 0.5
    title = r"$M_s$=62.6 GeV"
  if par=="lhs_ms":
    nx = 2
    ny = 0
    xlabel = r'$M_s$'
    ylabel = r'$\lambda_{HS}$'
    xmin = 10
    xmax = 120
    ymin = 0.1
    ymax = 0.5
    title = r"$\lambda_{S}$=0.1"
  if par=="ls_ms":
    nx = 2
    ny = 1
    xlabel = r'$M_s$'
    ylabel = r'$\lambda_{S}$'
    xmin = 10
    xmax = 120
    ymin = 0
    ymax = 0.2
    title = r"$\lambda_{HS}=0.3$"
  if par=="lhs_ls_ms":
    nx = 2
    ny = 1
    xlabel = r'$M_s$'
    ylabel = r'$\lambda_{S}$'
    xmin = 10
    xmax = 120
    ymin = 0
    ymax = 0.2
    title = r"$\lambda_{HS}\in[0.1,0.5]$"
    
  map1 = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,3], cmap=cm, vmin=vmin, vmax=vmax, alpha=0.8)
  if len(d1) > 0:
    ax.scatter(d1[:,nx], d1[:,ny], c='mediumseagreen', alpha=0.9, label=r'Only exist in BK scenario')
  if len(d2) > 0:
    ax.scatter(d2[:,nx], d2[:,ny], c='deepskyblue', alpha=0.9, label=r'Only exist in [' + compare_label+']')

  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
  ax.set_xlim(xmin,xmax)
  ax.set_ylim(ymin,ymax)
  ax.set_title(title)

  if par=="lhs_ls":
    ax.legend(loc=4)
    return map1

fig, axs = plt.subplots(2,2, figsize=(10, 10))

map1 = make_plot(axs[0,0], par="lhs_ls")
make_plot(axs[0,1], par="lhs_ms")
make_plot(axs[1,0], par="ls_ms")
make_plot(axs[1,1], par="lhs_ls_ms")

plt.subplots_adjust(hspace=0.4)

fig.suptitle(suptitle)


clb = plt.colorbar(map1, ax=axs, label="Relative error " + (r'$\Delta \gamma_{\rm EW} / \gamma_{\rm EW}^{bk}$' if show_gamma else r'$\Delta T_C/T_C^{BK}$' ))
#fig.tight_layout()

plt.savefig('compare_'+figure_name + ('_gamma' if show_gamma else '_Tc') + '.png')
plt.show()

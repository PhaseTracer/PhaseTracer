import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from plot_fun import fun_gamma, fun_diff, loaddata


#if len(sys.argv)<3 :
#  print "Use compare.py 0 1"
#  sys.exit()

#show_gamma = sys.argv[1] == '1'

#if sys.argv[2] == '1':
#  compare_label = r'$Q=2m_t$'
#  compare_name = "MSbar_ArnoldEspinosa_2mt_NoTreeEWSB"
#  figure_name = "2mt"
#elif sys.argv[2] == '2':
#  compare_label = r'Goldstone resummation'
#  compare_name = "MSbar_ArnoldEspinosa_mt_TreeEWSB"
#  figure_name = "Goldstone"
#elif sys.argv[2] == '3':
#  compare_label = r'Parwani resummation'
#  compare_name = "MSbar_Parwani_mt_NoTreeEWSB"
#  figure_name = "Parwani"
#elif sys.argv[2] == '4':
#  compare_label = r'$\xi=1$'
#  compare_name = "MSbar_ArnoldEspinosa_mt_NoTreeEWSB_xi"
#  

#label=r'On-shell like')
#label=r'High-$T$ approximation')
#label=r'PRM')
  
show_gamma = False
cm = 'autumn_r'
  
  
figure_name = "xi"

data1 = np.loadtxt("../2d_scan/lhs_ls_default.txt")
data2 = np.loadtxt("../2d_scan/lhs_ls_xi3.txt")
diff, d1, d2 = fun_diff(data2, data1, data1)
vmin = min(diff[:,3])
vmax = max(diff[:,3])

def make_plot(ax, par):
  #############################
  data1 = np.loadtxt("../2d_scan/"+par+"_default.txt")
  data2 = np.loadtxt("../2d_scan/"+par+"_xi3.txt")

  diff, d1, d2 = fun_diff(data2, data1, data1)
  
  # 0=lambda_hs, 1=lambda_s, 2=ms
  if par=="lhs_ls":    
    nx = 1
    ny = 2
    xlabel = r'$\lambda_{S}$'
    ylabel = r'$\lambda_{HS}$'
    xmin = 0
    xmax = 0.2
    ymin = 0.1
    ymax = 0.5
    title = r"$M_s$=62.6 GeV"
    
  if par=="ms_lhs":
    nx = 0
    ny = 2
    xlabel = r'$M_s$'
    ylabel = r'$\lambda_{HS}$'
    xmin = 10
    xmax = 120
    ymin = 0.1
    ymax = 0.5
    title = r"$\lambda_{S}$=0.1"
  if par=="ms_ls":
    nx = 0
    ny = 1
    xlabel = r'$M_s$'
    ylabel = r'$\lambda_{S}$'
    xmin = 10
    xmax = 120
    ymin = 0
    ymax = 0.2
    title = r"$\lambda_{HS}=0.3$"
    
  map1 = ax.scatter(diff[:,nx], diff[:,ny], c=diff[:,3], cmap=cm, vmin=vmin, vmax=vmax, alpha=0.8)

#  if len(d1) > 0:
#    ax.scatter(d1[:,nx], d1[:,ny], c='mediumseagreen', alpha=0.9, label=r'Only exist in BK scenario')
#  if len(d2) > 0:
#    ax.scatter(d2[:,nx], d2[:,ny], c='deepskyblue', alpha=0.9, label=r'Only exist in ')

  ax.set_xlabel(xlabel)
  ax.set_ylabel(ylabel)
#  ax.set_xlim(xmin,xmax)
#  ax.set_ylim(ymin,ymax)
  ax.set_title(title)

  clb = plt.colorbar(map1, ax=axs, label=r"$T_C^{(\xi=3)}-T_C^{(\xi=0)}$" )

#  if par=="lhs_ls":
##    ax.legend(loc=4)
#    return map1

fig, axs = plt.subplots(1,3, figsize=(15, 5))

make_plot(axs[0], par="lhs_ls")
make_plot(axs[1], par="ms_lhs")
make_plot(axs[2], par="ms_ls")

#fig.tight_layout()

#plt.subplots_adjust(hspace=0.2)



plt.savefig('compare_'+figure_name + ('_gamma' if show_gamma else '_Tc') + '.png')
plt.show()

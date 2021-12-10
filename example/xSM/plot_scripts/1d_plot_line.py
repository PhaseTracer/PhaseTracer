import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
#sys.path.append( '../' )
from plot_fun import fun_gamma, fun_diff, loaddata

plot_methods = False
plot_Goldstone = False
plot_xi = False
plot_xi_zoomin = False 

plot_scale = True

figure_format = "pdf"

cmap = cm.get_cmap('rainbow')

def plot_for_1d(data, x_num, label, column):

  sel = data[:,3]>0

  x = data[:,x_num]
  TC = data[:,4]
  gamma = fun_gamma(data)

  ax = axs[0,column]
  ax.plot(x[sel], TC[sel], label=label, alpha=1)

  ax = axs[1,column]
  ax.plot(x[sel], gamma[sel], label=label, alpha=1)

if plot_methods or plot_Goldstone:
  if plot_methods:
    
#    names = [ ["default", "default"],
#              ["2mt", r"$Q=2m_t$"],
#              ["OSlike", r"OS-like"],
#              ["xi3", r"$\xi=3$"],
#              ["HT", r"HT"],
#              ["PRM","PRM"],
#              ["covariant_gauge", r"covariant, $\xi=3$"]]

#    names = [ ["default", "MSbar, ArnoldEspinosa"],
#              ["OSlike", r"OS-like, ArnoldEspinosa"],
#              ["Parwani","MSbar, Parwani"],
#              ["noDaisy", r"MSbar, no daisy"]] 
              
    names = [ ["default", "MSbar"],
              ["OSlike", r"OS-like"],
              ["PRM_woFS_0L", r"PRM(0-L)"],
              ["PRM_woFS", r"PRM(1-L)"],
              ["HT","HT"]] 
              
    ncolumn = 3
  if plot_Goldstone:
    names = [ ["default", "Goldstone resummation"],
              ["1L_EWSB", "One-loop EWSB"] ]
    ncolumn = 4

  fig, axs = plt.subplots(2, ncolumn, figsize=(13, 6))

  for name in names:
    plot_for_1d(np.loadtxt("../1d_bks/lambda_hs_"+name[0]+".txt"), 2, name[1], 0)
    plot_for_1d(np.loadtxt("../1d_bks/lambda_s_"+name[0]+".txt"), 1, name[1], 1)
    plot_for_1d(np.loadtxt("../1d_bks/m_s_"+name[0]+".txt"), 0, name[1], 2)
  
  if plot_Goldstone:
    plot_for_1d(np.loadtxt("../1d_bks/Rxi_MSbar_resummation.txt"), 10, "Goldstone resummation", 3)
    plot_for_1d(np.loadtxt("../1d_bks/Rxi_MSbar_1L_EWSB.txt"), 10, "One-loop EWSB", 3)
    plot_for_1d(np.loadtxt("../1d_bks/Rxi_MSbar_no.txt"), 10, "Nothing", 3)

if plot_scale:
  names = [ ["PRM_0L_noRGE_mt", "_PRM_0L_noRGE_mt"],
            ["PRM_woFS_noRGE_mt", "PRM_woFS_noRGE_mt"]
               ]
  ncolumn = 3
  fig, axs = plt.subplots(2, ncolumn, figsize=(13, 6))
  for name in names:
    plot_for_1d(np.loadtxt("../1d_bks/lambda_hs_"+name[0]+".txt"), 2, name[1], 0)
    plot_for_1d(np.loadtxt("../1d_bks/lambda_s_"+name[0]+".txt"), 1, name[1], 1)
    plot_for_1d(np.loadtxt("../1d_bks/m_s_"+name[0]+".txt"), 0, name[1], 2)
    

if plot_xi:
  
  fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    
  names = [ ["MSbar", r"$\overline{\rm MS}$(Resum $M_G^2$)"],
            ["MSbar_1L_EWSB", r"$\overline{\rm MS}$(Modified $\mu$)"],
            ["MSbar_no", r"$\overline{\rm MS}$(None)"],
            ["PRM", "PRM(1L tadpole)"],
            ["PRM_0L", r"PRM"]]
  ncolumn = 2
    
  for name in names:
    plot_for_1d(np.loadtxt("../1d_bks/Rxi_"+name[0]+".txt"), 10, name[1], 0)
    plot_for_1d(np.loadtxt("../1d_bks/covariant_"+name[0]+".txt"), 10, name[1], 1)

#  plot_for_1d(np.loadtxt("../1d_bks/Rxi_PRM_FS_0L.txt"), 10, "PRM(FS,0-L)", 0)
  
  data = np.loadtxt("../1d_bks/Rxi_HT.txt")
  TC = data[:,4]
  gamma = fun_gamma(data)
  axs[0,0].plot([0,10], [TC[0],TC[-1]], label="HT", alpha=1)
  axs[0,1].plot([0,10], [TC[0],TC[-1]], label="HT", alpha=1)
  axs[1,0].plot([0,10], [gamma[0],gamma[-1]], label="HT", alpha=1)
  axs[1,1].plot([0,10], [gamma[0],gamma[-1]], label="HT", alpha=1)

  axs[0,0].set_title(r"R-$\xi$ gauge")
  axs[0,1].set_title(r"Covariant gauge")



  
for ii in range(2):
    for jj in range(ncolumn):
      axs[ii,jj].grid(axis='x', alpha=0.75)
      axs[ii,jj].grid(axis='y', alpha=0.75)
      
      if ii == 0:
        axs[ii,jj].set_ylabel(r"$T_C$ (GeV)")
#        if jj<1:
#          axs[ii,jj].legend(loc=3)
#        else:
#          axs[ii,jj].legend(loc=4)

      else:
        axs[ii,jj].set_ylim(0,5)
        axs[ii,jj].set_ylabel(r"$\gamma_{\rm EW}$")
#        if jj<1:
#          axs[ii,jj].legend(loc=2)
#        else:
#          axs[ii,jj].legend(loc=1)
      
      if plot_methods or plot_Goldstone or plot_scale:
        if jj == 0:
          axs[ii,jj].set_xlabel(r"$\lambda_{hs}$")
        elif jj == 1:
          axs[ii,jj].set_xlabel(r"$\lambda_{s}$")
        elif jj == 2:
          axs[ii,jj].set_xlabel(r"$m_{s}$ (GeV)")
        else:
          axs[ii,jj].set_xlabel(r"$\xi$")
      else:
        if jj == 0:
          axs[ii,jj].set_xlabel(r"$\xi$")
        else:
          axs[ii,jj].set_xlabel(r"$\xi_W=\xi_B$")
        if plot_xi:
          axs[ii,jj].set_xlim(0,10)
        if ii == 1:
          if plot_xi:
            axs[ii,jj].set_ylim(1.5,2.5)
        if plot_xi_zoomin:
          axs[ii,jj].set_xlim(0,0.5)
          if ii == 1:
            axs[ii,jj].set_ylim(1.75,2.05)
          else:
            axs[ii,jj].set_ylim(110,114)

if plot_xi:
  axs[0,1].legend(bbox_to_anchor=(1,0), loc=3)
if plot_scale:
  axs[0,0].legend(loc=3)

fig.tight_layout()

if plot_methods:
  plt.savefig('1d_methods.png')
elif plot_Goldstone:
  plt.savefig('1d_Goldstone.png')
elif plot_xi:
  if not plot_xi_zoomin:
    plt.savefig('1d_xi.'+figure_format)
  else:
    plt.savefig('1d_xi_zoom_in.'+figure_format)
plt.show()

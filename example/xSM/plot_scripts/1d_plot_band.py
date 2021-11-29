import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
#sys.path.append( '../' )
from plot_fun import fun_gamma, fun_diff, loaddata

cmap = cm.get_cmap('rainbow')

plot_scale = True
plot_xi = False

def line_for_1d(data, x_num, color, column, label = ""):
  sel = data[:,3]>0

  x = data[:,x_num][sel]
  TC = data[:,4][sel]
  gamma = fun_gamma(data)[sel]

  ax = axs[0,column]
  if label == "":
    ax.plot(x, TC, color=color, alpha=1)
  else:
    ax.plot(x, TC, color=color, alpha=1, label=label)
    
  ax = axs[1,column]
  if label == "":
    ax.plot(x, gamma, color=color, alpha=1)
  else:
    ax.plot(x, gamma, color=color, alpha=1, label=label)
    
    
def range_for_1d(data1, data2, x_num, label, column, color):

  sel1 = data1[:,3]>0
  x1 = data1[:,x_num][sel1]
  TC1 = data1[:,4][sel1]
  gamma1 = fun_gamma(data1)[sel1]

  sel2 = data2[:,3]>0
  x2 = data2[:,x_num][sel2]
  TC2 = data2[:,4][sel2]
  gamma2 = fun_gamma(data2)[sel2]

  fTC1 = interp1d(x1, TC1, kind='linear')
  fTC2 = interp1d(x2, TC2, kind='linear')
  fgamma1 = interp1d(x1, gamma1, kind='linear')
  fgamma2 = interp1d(x2, gamma2, kind='linear')
  
  x = np.linspace(max(min(x1),min(x2)), min(max(x1),max(x2)), num=50, endpoint=True)

  ax = axs[0,column]
  ax.fill_between(x,fTC1(x),fTC2(x), color=color, alpha=0.3, linewidth=0, label=label)
  
  ax = axs[1,column]
  ax.fill_between(x,fgamma1(x),fgamma2(x), color=color, alpha=0.3, linewidth=0, label=label)

fig, axs = plt.subplots(2, 3, figsize=(12, 6))

for name in [["m_s",0,2], ["lambda_s",1,1], ["lambda_hs",2,0],]:


  if plot_scale:

#    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_noRGE_mt.txt"), name[1], "b", name[2])
#    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_noRGE_05mt.txt"),
#                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_noRGE_2mt.txt"),
#                 name[1], r"PRM(no RGE), $Q\in[m_t/2,2m_t]$", name[2], 'b')

#    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noRGE_mt.txt"), name[1], "purple", name[2])
#    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noRGE_05mt.txt"),
#                 np.loadtxt("../1d_bks/"+name[0]+"_noRGE_2mt.txt"),
#                 name[1], r"MS(no RGE), $Q\in[m_t/2,2m_t]$", name[2], 'purple')

    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_mt.txt"), name[1], "k", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_2mt.txt"),
                 name[1], r"PRM(No RGE), $Q\in[m_t/2,2m_t]$", name[2], 'k')
                 
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_mt.txt"), name[1], "yellow", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_2mt.txt"),
                 name[1], r"PRM(No RGE, FS), $Q\in[m_t/2,2m_t]$", name[2], 'yellow')

    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_mt.txt"), name[1], "r", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_2mt.txt"),
                 name[1], r"PRM, $Q\in[m_t/2,2m_t]$", name[2], 'r')
                 
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_mt.txt"), name[1], "g", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_2mt.txt"),
                 name[1], r"MS, $Q\in[m_t/2,2m_t]$", name[2], 'g')
                 
                 
                 
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noD_mt.txt"), name[1], "b", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noD_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_noD_2mt.txt"),
                 name[1], r"MS(no daisy), $Q\in[m_t/2,2m_t]$", name[2], 'b')




                 
  elif plot_xi:     
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_xi1.txt"), name[1], "g", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_xi0.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_xi3.txt"),
                 name[1], r"MS, $xi\in[0,3]$", name[2], 'g')
                 
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noD_xi1.txt"), name[1], "r", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noD_xi0.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_noD_xi3.txt"),
                 name[1], r"MS(no daisy), $xi\in[0,3]$", name[2], 'r')

    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_xi1.txt"), name[1], "b", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_xi0.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_xi3.txt"),
                 name[1], r"PRM(1L), $xi\in[0,3]$", name[2], 'b')

    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_0L.txt"), name[1], "k", name[2], r"PRM(0L)")
               
  
for ii in range(2):
  for jj in range(3):
    axs[ii,jj].grid(axis='x', alpha=0.75)
    axs[ii,jj].grid(axis='y', alpha=0.75)
    
    if ii == 0:
      axs[ii,jj].set_ylabel(r"$T_C$ (GeV)")
#      if jj<1:
#        axs[ii,jj].legend(loc=3)
#      else:
#        axs[ii,jj].legend(loc=4)
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

if plot_scale:
  plt.savefig('1d_scale.png')
elif plot_xi:
  plt.savefig('1d_xi.png')

plt.show()

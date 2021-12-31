import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.interpolate import interp1d
#sys.path.append( '../' )
from plot_fun import fun_gamma, fun_diff, loaddata

cmap = cm.get_cmap('rainbow')

plot_scale = True
plot_xi = False

def selection(for_TC, x_num, data):
  if for_TC:
    sel = (data[:,3]>0) & (data[:,4]<200.)
    y = data[:,4][sel]
  else:
    sel = (data[:,3]>0) & (fun_gamma(data)>0)
    y = fun_gamma(data)[sel]
  x = data[:,x_num][sel]
  return x,y

def line_for_1d(data, x_num, color, column, label = "", linestyle="-", alpha=1):
  x_TC, TC = selection(True, x_num, data)
  x_gamma, gamma = selection(False, x_num, data)
  
  ax = axs[0,column]
  if label == "":
    ax.plot(x_TC, TC, color=color, linestyle=linestyle, alpha=alpha)
  else:
    ax.plot(x_TC, TC, color=color, linestyle=linestyle, alpha=alpha, label=label)
    
  ax = axs[2,column]
  if label == "":
    ax.plot(x_gamma, gamma, color=color, linestyle=linestyle, alpha=alpha)
  else:
    ax.plot(x_gamma, gamma, color=color, linestyle=linestyle, alpha=alpha, label=label)
    
def interpolate(x1, y1, x2, y2):
  if min(x1) < min(x2):
    x2 = np.append(x1[0],x2)
    y2 = np.append(y1[0],y2)
  elif min(x1) > min(x2):
    x1 = np.append(x2[0],x1)
    y1 = np.append(y2[0],y1)
  if max(x1) > max(x2):
    x2 = np.append(x2,x1[-1])
    y2 = np.append(y2,y1[-1])
  elif max(x1) < max(x2):
    x1 = np.append(x1,x2[-1])
    y1 = np.append(y1,y2[-1])

  f1 = interp1d(x1, y1, kind='linear')
  f2 = interp1d(x2, y2, kind='linear')
  return f1, f2  
  
def range_for_1d(data1, data2, x_num, label, column, color):

  x_TC1, TC1 = selection(True, x_num, data1)
  x_TC2, TC2 = selection(True, x_num, data2)

  x_gamma1, gamma1 = selection(False, x_num, data1)
  x_gamma2, gamma2 = selection(False, x_num, data2)

  fTC1, fTC2 = interpolate(x_TC1, TC1, x_TC2, TC2)
  fgamma1, fgamma2 = interpolate(x_gamma1, gamma1, x_gamma2, gamma2)
  
  x_TC = np.linspace(min(min(x_TC1),min(x_TC2)), max(max(x_TC1),max(x_TC2)), num=100, endpoint=True)
  x_gamma = np.linspace(min(min(x_gamma1),min(x_gamma2)), max(max(x_gamma1),max(x_gamma2)), num=100, endpoint=True)

  ax = axs[0,column]
  ax.fill_between(x_TC,fTC1(x_TC),fTC2(x_TC), color=color, alpha=0.3, linewidth=0, label=label)

  ax = axs[2,column]
  ax.fill_between(x_gamma,fgamma1(x_gamma),fgamma2(x_gamma), color=color, alpha=0.3, linewidth=0, label=label)
  
#  diff_TC = fun_diff(data1, data2, data1, use_abs=True, sort=False)
#  ax = axs[1,column]
#  ax.plot(diff_TC[:,x_num], diff_TC[:,4], color=color, label=label)
  x_TC = np.linspace(max(min(x_TC1),min(x_TC2)), min(max(x_TC1),max(x_TC2)), num=100, endpoint=True)
  ax = axs[1,column]
  ax.plot(x_TC,(fTC1(x_TC)-fTC2(x_TC)), color=color, label=label)
  
  
#  diff_gamma = fun_diff(data1, data2, data1, show_gamma=True, use_abs=True, sort=False)
#  ax = axs[3,column]
#  ax.plot(diff_gamma[:,x_num], diff_gamma[:,4], color=color, label=label)
  x_gamma = np.linspace(max(min(x_gamma1),min(x_gamma2)), min(max(x_gamma1),max(x_gamma2)), num=100, endpoint=True)
  ax = axs[3,column]
  ax.plot(x_gamma,(fgamma1(x_gamma)-fgamma2(x_gamma)), color=color, label=label)

  
fig, axs = plt.subplots(4, 3, figsize=(12, 12))

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


    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_mt.txt"), name[1], "g", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_2mt.txt"),
                 name[1], r"$\overline{\rm MS}$", name[2], 'g')
                 
                
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noD_mt.txt"), name[1], "b", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noD_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_noD_2mt.txt"),
                 name[1], r"$\overline{\rm MS}$(no daisy)", name[2], 'b')

    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noRGE_mt.txt"), name[1], "gray", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_noRGE_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_noRGE_2mt.txt"),
                 name[1], r"$\overline{\rm MS}$(no RGE)", name[2], 'gray')

    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_mt.txt"), name[1], "r", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_2mt.txt"),
                 name[1], r"PRM", name[2], 'r')
                                  
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_mt.txt"), name[1], "purple", name[2])
    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_05mt.txt"),
                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_2mt.txt"),
                 name[1], r"PRM(no RGE)", name[2], 'purple')
    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_0L_noRGE_2mt.txt"), name[1], "purple", name[2], linestyle= "--", alpha=0.4)
                 
#    line_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_mt.txt"), name[1], "k", name[2])
#    range_for_1d(np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_05mt.txt"),
#                 np.loadtxt("../1d_bks/"+name[0]+"_PRM_woFS_noRGE_2mt.txt"),
#                 name[1], r"PRM(No RGE woFS), $Q\in[m_t/2,2m_t]$", name[2], 'k')

                 


                 
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
               
  
for ii in range(4):
  for jj in range(3):
    axs[ii,jj].grid(axis='x', alpha=0.5)
    axs[ii,jj].grid(axis='y', alpha=0.5)
    
    if ii == 0:
      axs[ii,jj].set_ylim(40,170)
      axs[ii,jj].set_ylabel(r"$T_C$ (GeV)")
    elif ii == 2:
      axs[ii,jj].set_ylim(0,6)
      axs[ii,jj].set_ylabel(r"$\gamma_{\rm EW}$")
    elif ii == 1:
      axs[ii,jj].set_ylim(-20,30)
      axs[ii,jj].set_ylabel(r"$T_C^{(Q=\frac{1}{2}m_t)}-T_C^{(Q=2m_t)}$ (GeV)")
    elif ii == 3:
      axs[ii,jj].set_ylim(-4,2)
      axs[ii,jj].set_ylabel(r"$\gamma^{(Q=\frac{1}{2}m_t)}-\gamma^{(Q=2m_t)}$ (GeV)")
      
    if jj == 0:
      axs[ii,jj].set_xlabel(r"$\lambda_{hs}$")
      axs[ii,jj].set_xlim(0.1,0.4)
    elif jj == 1:
      axs[ii,jj].set_xlabel(r"$\lambda_{s}$")
      axs[ii,jj].set_xlim(0.05,0.2)
    else:
      axs[ii,jj].set_xlabel(r"$m_{s}$ (GeV)")
      axs[ii,jj].set_xlim(40,100)

    axs[3,0].legend(loc=2)
    axs[2,0].legend(loc=2)

fig.tight_layout()

if plot_scale:
  plt.savefig('1d_scale.pdf')
elif plot_xi:
  plt.savefig('1d_xi.png')

plt.show()

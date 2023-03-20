#!/usr/bin/env python
"""
Plot phases
===========
"""

import re
import os
import sys
from scipy.interpolate import interp1d
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
from matplotlib.legend_handler import HandlerPatch
from style import style


ARROW_DELTA_T = 10.
QUIVER_ARROW = {"zorder": 10, "angles": 'xy', "scale_units": 'xy', "scale": 1, "units": 'dots', "width": 4, "minlength": 3.}

GEV = ""
FONTSIZE = "x-small"


LINE = {"lw": 2, "alpha": 1}
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

style()
fig = plt.figure(figsize=(12, 5))
name = ["86.500000", "173.000000", "346.000000"]
label = [r"\frac12 m_t", "m_t", "2 m_t"]

def set_label(jj, T, title=""):
  if jj == 0:
    plt.xlabel(r"$\phi_h$ [GeV]")
    t = r"$\phi_s=0$"
  else:
    plt.xlabel(r"$\phi_s$ [GeV]")
    t = r"$\phi_h=0$"
  plt.ylabel(r"$V_{\rm eff}$ [GeV${}^4$]")
  plt.title(r"$T="+str(T)+"$, "+t+", "+title)

def find_loc(data, value):
  diff = abs(data-value)
  return diff == min(diff)


for ii in range(3):
  T = 0
  
  tree_min = [242.9644987251219, 259.672598075815, 276.555894966469]  
  data = np.loadtxt("potential_line/"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(1,2,1)

  
  min_loc = find_loc(data[:,2],tree_min[ii])
  plt.scatter(data[min_loc][:,2], data[min_loc][:,1], marker="+", s=80, c="black",
    label="Tree-level minimum" if ii == 0 else None, zorder=10)
  
  VEV = data[:,1] == min(data[:,1])
  plt.scatter(data[VEV][:,2], data[VEV][:,1], marker="x", s=40, c="black",
      label="One-loop minimum" if ii == 0 else None, zorder=10)

  plt.plot(data[:,2], data[:,1], label = "$Q="+label[ii]+"$", c=colors[ii], **LINE)
    
  set_label(0, T, "w RGE")
  plt.legend(loc=2)
  plt.xlim(237,282)
  plt.ylim(-1.3e8, -1.15e8)

  tree_min = [243.3873771875804, 259.3881307737823, 278.9974402653135]
  data = np.loadtxt("potential_line/no_RGE_"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(1,2,2)
  
  plt.plot(data[:,2], data[:,1], label = "$Q="+name[ii][0:4]+"$ GeV", c=colors[ii], **LINE)
  
  min_loc = find_loc(data[:,2],tree_min[ii])
  plt.scatter(data[min_loc][:,2], data[min_loc][:,1], marker="+", s=80, c="black", zorder=10)
  
  VEV = data[:,1] == min(data[:,1])
  plt.scatter(data[VEV][:,2], data[VEV][:,1], marker="x", s=40, c="black", zorder=10)
  
  set_label(0, T, "w/o RGE")

  plt.xlim(237,282)
  plt.ylim(-1.3e8, -1.15e8)
  
 

fig.tight_layout()

fig.subplots_adjust(wspace=0.25,hspace=0.5)

plt.savefig("potential_scale.pdf")
plt.show()











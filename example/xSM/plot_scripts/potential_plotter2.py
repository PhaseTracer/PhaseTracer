#!/usr/bin/env python
"""
Plot phases
===========
"""

import re
import os
import sys
from StringIO import StringIO
from scipy.interpolate import interp1d
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
from matplotlib.legend_handler import HandlerPatch

ARROW_DELTA_T = 10.
QUIVER_ARROW = {"zorder": 10, "angles": 'xy', "scale_units": 'xy', "scale": 1, "units": 'dots', "width": 4, "minlength": 3.}

GEV = ""
FONTSIZE = "x-small"


LINE = {"lw": 1, "alpha": 0.8}
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fig = plt.figure(figsize=(9,5))
name = ["86.500000", "173.000000", "346.000000"]

def set_label(jj, T, title=""):
  if jj == 0:
    plt.xlabel(r"$\phi_h$ [GeV]")
    t = r"$\phi_s=0$"
  else:
    plt.xlabel(r"$\phi_s$ [GeV]")
    t = r"$\phi_h=0$"
  plt.ylabel(r"$V_{\rm eff}$")
  plt.title(r"T="+str(T)+" GeV, "+t+", "+title)

def find_loc(data, value):
  diff = abs(data-value)
  return diff == min(diff)


for ii in range(3):
  T = 0
  
  tree_min = [242.98498, 259.62229, 276.41735]  
  data = np.loadtxt("potential_line/"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(1,2,1)
  plt.plot(data[:,2], data[:,1], label = "Q="+name[ii][0:4]+" GeV", c=colors[ii], **LINE)
  
  min_loc = find_loc(data[:,2],tree_min[ii])
  plt.plot(data[min_loc][:,2], data[min_loc][:,1], marker="+", markersize=20, c=colors[ii])
  
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,2], data[VEV][:,1], marker="x", markersize=10, c=colors[ii])
  
  set_label(0, T, "w RGE")
  plt.legend(loc=2)
  plt.xlim(237,282)
  
  tree_min = [243.3731362896345, 259.4444281087366, 279.1590955719373]  
  data = np.loadtxt("potential_line/no_RGE_"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(1,2,2)
  
  plt.plot(data[:,2], data[:,1], label = "Q="+name[ii][0:4]+" GeV", c=colors[ii], **LINE)
  
  min_loc = find_loc(data[:,2],tree_min[ii])
  plt.plot(data[min_loc][:,2], data[min_loc][:,1], marker="+", markersize=20, c=colors[ii])
  
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,2], data[VEV][:,1], marker="x", markersize=10, c=colors[ii])
  
  set_label(0, T, "w/o RGE")

  plt.xlim(237,282)
  
 

fig.tight_layout()

fig.subplots_adjust(wspace=0.25,hspace=0.5)

plt.savefig("potential_scale.pdf")
plt.show()











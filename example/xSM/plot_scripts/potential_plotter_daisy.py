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

ARROW_DELTA_T = 10.
QUIVER_ARROW = {"zorder": 10, "angles": 'xy', "scale_units": 'xy', "scale": 1, "units": 'dots', "width": 4, "minlength": 3.}

GEV = ""
FONTSIZE = "x-small"


LINE = {"lw": 1, "alpha": 0.8}
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fig = plt.figure(figsize=(9,10))
name = ["Arnold Espinosa", "Parwani", "no daisy"]

def set_label(jj, T):
  if jj == 0:
    plt.xlabel(r"$\phi_h$ [GeV]")
    t = r"$\phi_s=0$"
  else:
    plt.xlabel(r"$\phi_s$ [GeV]")
    t = r"$\phi_h=0$"
  plt.ylabel(r"$V_{\rm eff}/T^4$")
  plt.title(r"T="+str(T)+" GeV, "+t)

for ii in range(3):
  T = 160
  data = np.loadtxt("potential_line/"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(3,2,1)
  plt.plot(data[:,2], data[:,1]/pow(T,4), label = name[ii], c=colors[ii], **LINE)
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,2], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
  set_label(0, T)
  plt.legend(loc=2)
  plt.xlim(-50,300)
  data = np.loadtxt("potential_line/"+str(T)+"_2_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(3,2,2)
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,3], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
  set_label(1, T)
  plt.xlim(-50,300)

  T = 110
  data = np.loadtxt("potential_line/"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(3,2,3)
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,2], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
  set_label(0, T)
  plt.xlim(-50,300)
  data = np.loadtxt("potential_line/"+str(T)+"_2_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(3,2,4)
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,3], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
  set_label(1, T)
  plt.xlim(-50,300)

  T = 10
  data = np.loadtxt("potential_line/"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(6,4,17)
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  set_label(0, T)
  plt.xlim(-0.5,2.5)
  data = np.loadtxt("potential_line/"+str(T)+"_2_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(6,4,22)
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,2], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
  set_label(0, T)
  plt.xlim(244.5,250.5)
  data = np.loadtxt("potential_line/"+str(T)+"_3_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(6,4,19)
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  set_label(1, T)
  plt.xlim(-0.5,2.5)
  data = np.loadtxt("potential_line/"+str(T)+"_4_"+name[ii].replace(" ", "")+"_potential_line.dat")
  plt.subplot(6,4,24)
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  VEV = data[:,1] == min(data[:,1])
  plt.plot(data[VEV][:,3], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
  set_label(1, T)
  plt.xlim(217.5,223.5)

fig.tight_layout()

fig.subplots_adjust(wspace=0.25,hspace=0.5)

plt.savefig("potential_daisy.pdf")
plt.show()











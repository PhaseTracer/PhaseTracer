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


#folder = "../../../build/"
folder = "potential_line/"

fig = plt.figure(figsize=(9,10))
name = ["Arnold Espinosa", "Parwani", "no daisy"]


m_s = 65
lambda_hs = 0.3
lambda_s = 0.1
v = 247.4534576428087
g = 0.6477108911331643
gp = 0.3585642676423744
yt = 0.9341420067357675
yb = 0.01547368604435034
ytau = 0.01001414937355005
mh = 125

def square(x):
  return x*x
def pow_4(x):
  return x*x*x*x
    
muh_sq = -0.5 * mh * mh
lambda_h = -muh_sq / square(v)
mus_sq = square(m_s) - lambda_hs * square(v) / 2.
g_sq = g*g
gp_sq = gp*gp
yt_sq = yt*yt

def HT(phi,T):
  muh_sq_T = muh_sq + square(T) / 48. * \
        (9. * g_sq + 3. * gp_sq + 12. * yt_sq + 24. * lambda_h + 2. * lambda_hs)
  mus_sq_T = mus_sq + square(T) / 12. * (2. * lambda_hs + 3. * lambda_s)
    
  return 0.5 * muh_sq_T * square(phi[0]) + \
           0.25 * lambda_h * pow_4(phi[0]) + \
           0.25 * lambda_hs * square(phi[0]) * square(phi[1]) + \
           0.5 * mus_sq_T * square(phi[1]) + \
           0.25 * lambda_s * pow_4(phi[1])

def HT_h(x,T):
  phi = [x, 0]
  return HT(phi,T)

def HT_s(x,T):
  phi = [0, x]
  return HT(phi,T)

def set_label(jj, T):
  if jj == 0:
    plt.xlabel(r"$\phi_h$ [GeV]")
    t = r"$\phi_s=0$"
  else:
    plt.xlabel(r"$\phi_s$ [GeV]")
    t = r"$\phi_h=0$"
  if T==0:
    plt.ylabel(r"$V_{\rm eff}$")
  else:
    plt.ylabel(r"$V_{\rm eff}/T^4$")
  plt.title(r"T="+str(T)+" GeV, "+t)

for ii in range(3):
#  T = 150
#  data = np.loadtxt("potential_line/"+str(T)+"_1_"+name[ii].replace(" ", "")+"_potential_line.dat")
#  plt.subplot(3,2,1)
#  plt.plot(data[:,2], data[:,1]/pow(T,4), label = name[ii], c=colors[ii], **LINE)
#  VEV = data[:,1] == min(data[:,1])
#  plt.plot(data[VEV][:,2], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
#  set_label(0, T)
#  plt.legend(loc=2)
#  plt.xlim(-50,300)
#  data = np.loadtxt("potential_line/"+str(T)+"_2_"+name[ii].replace(" ", "")+"_potential_line.dat")
#  plt.subplot(3,2,2)
#  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
#  VEV = data[:,1] == min(data[:,1])
#  plt.plot(data[VEV][:,3], data[VEV][:,1]/pow(T,4), ".", c=colors[ii])
#  set_label(1, T)
#  plt.xlim(-50,300)

  T = 150

  plt.subplot(3,2,1)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_Q_173.000000_potential_line.dat")
  if ii==0:
    plt.plot(data[:,2], HT_h(data[:,2], T)/pow(T,4), label = 'HT', c='k', linestyle="--", **LINE)
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], label = name[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_Q_346.000000_potential_line.dat")
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_Q_86.500000_potential_line.dat")
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  set_label(0, T)
  plt.xlim(-50,300)
  plt.legend(loc=2)

  plt.subplot(3,2,2)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_Q_173.000000_potential_line.dat")
  if ii==0:
    plt.plot(data[:,3], HT_s(data[:,3], T)/pow(T,4), label = 'HT', c='k', linestyle="--", **LINE)
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_Q_346.000000_potential_line.dat")
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_Q_86.500000_potential_line.dat")
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  set_label(1, T)
  plt.xlim(-50,300)


  T = 100

  plt.subplot(3,2,3)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_Q_173.000000_potential_line.dat")
  if ii==0:
    plt.plot(data[:,2], HT_h(data[:,2], T)/pow(T,4), label = 'HT', c='k', linestyle="--", **LINE)
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_Q_346.000000_potential_line.dat")
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_Q_86.500000_potential_line.dat")
  plt.plot(data[:,2], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  set_label(0, T)
  plt.xlim(-50,300)


  plt.subplot(3,2,4)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_Q_173.000000_potential_line.dat")
  if ii==0:
    plt.plot(data[:,3], HT_s(data[:,3], T)/pow(T,4), label = 'HT', c='k', linestyle="--", **LINE)
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_Q_346.000000_potential_line.dat")
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_Q_86.500000_potential_line.dat")
  plt.plot(data[:,3], data[:,1]/pow(T,4), c=colors[ii], linestyle="--", **LINE)
  set_label(1, T)
  plt.xlim(-50,300)



  T = 0
  plt.subplot(6,4,17)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_1_Q_173.000000_potential_line.dat")
  plt.plot(data[:,2], data[:,1], c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_1_Q_346.000000_potential_line.dat")
  plt.plot(data[:,2], data[:,1], c=colors[ii], linestyle="--", **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_1_Q_86.500000_potential_line.dat")
  plt.plot(data[:,2], data[:,1], c=colors[ii], linestyle="--", **LINE)
  set_label(0, T)
  plt.xlim(-0.5,10)
  plt.subplot(6,4,22)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_2_Q_173.000000_potential_line.dat")
  plt.plot(data[:,2], data[:,1], c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_2_Q_346.000000_potential_line.dat")
  plt.plot(data[:,2], data[:,1], c=colors[ii], linestyle="--", **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_1_2_Q_86.500000_potential_line.dat")
  plt.plot(data[:,2], data[:,1], c=colors[ii], linestyle="--", **LINE)
  set_label(0, T)
  plt.xlim(240.5,260.5)


  plt.subplot(6,4,19)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_1_Q_173.000000_potential_line.dat")
  plt.plot(data[:,3], data[:,1], c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_1_Q_346.000000_potential_line.dat")
  plt.plot(data[:,3], data[:,1], linestyle="--", c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_1_Q_86.500000_potential_line.dat")
  plt.plot(data[:,3], data[:,1], linestyle="--", c=colors[ii], **LINE)
  set_label(1, T)
  plt.xlim(-0.5,10)
  plt.subplot(6,4,24)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_2_Q_173.000000_potential_line.dat")
  plt.plot(data[:,3], data[:,1], c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_2_Q_346.000000_potential_line.dat")
  plt.plot(data[:,3], data[:,1], linestyle="--", c=colors[ii], **LINE)
  data = np.loadtxt(folder+name[ii].replace(" ", "")+"_"+str(T)+"_2_2_Q_86.500000_potential_line.dat")
  plt.plot(data[:,3], data[:,1], linestyle="--", c=colors[ii], **LINE)
  set_label(1, T)
  plt.xlim(215,230)

fig.tight_layout()

fig.subplots_adjust(wspace=0.25,hspace=0.5)

plt.savefig("potential_scale_daisy_2.pdf")
plt.show()











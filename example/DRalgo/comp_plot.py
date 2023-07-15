import os
import random
import numpy as np
import multiprocessing as multi
from scipy import integrate
import scipy.optimize as opt
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from cosmoTransitions.tunneling1D import SingleFieldInstanton


DR1 =np.loadtxt("1d_bks/TcList1.dat")
DR2 =np.loadtxt("1d_bks/TcList2.dat")

PT =np.loadtxt("1d_bks/lam.txt")


fig, axes = plt.subplots(1, 1, figsize=(12, 7))


axes.plot(PT[:,2], PT[:,4],c='r')
axes.plot(DR1[:,0], DR1[:,1])
axes.plot(DR2[:,0], DR2[:,1])

plt.tight_layout()
plt.show()

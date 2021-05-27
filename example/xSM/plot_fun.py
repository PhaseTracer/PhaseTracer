import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def fun_gamma(data):
  return abs(data[:,5]-data[:,7])/data[:,4]

def fun_gamma_line(line):
  return abs(line[5]-line[7])/line[4]
  
def loaddata(name, gamma_min=0.1):
  data = np.loadtxt(name)
  data = data[data[:,3]>0]
  return data[fun_gamma(data)>gamma_min]

def fun_diff(data1, data2, data0, show_gamma=False, gamma_min = 0.7):
  if len(data1) != len(data2) or len(data1) != len(data0):
    print "Wrong data files"
    sys.exit()
  data_diff = []
  for ii in range(len(data1)):
    if data1[ii][0] != data2[ii][0] or data1[ii][0] != data0[ii][0]:
      print "Wrong data files"
      sys.exit()
    ms = data1[ii][0]
    lambda_s = data1[ii][1]
    lambda_hs = data1[ii][2]

    size_1 = data1[ii][3]
    size_2 = data2[ii][3]
    size_0 = data0[ii][3]
    
    gamma_1 = fun_gamma_line(data1[ii])
    gamma_2 = fun_gamma_line(data2[ii])
    gamma_0 = fun_gamma_line(data0[ii])

    if size_1 > 0 and size_2 > 0 and size_0>0:
      if gamma_0 > gamma_min and gamma_0 < 5: 
        if show_gamma:
          diff = abs(gamma_1 - gamma_2)
          central = gamma_0
        else:
          diff = abs(data1[ii][4] - data2[ii][4])
          central = data0[ii][4]
        data_diff.append([ms, lambda_s, lambda_hs, diff/central, central])

  data_diff.sort(key=(lambda x:x[3])) 
  
  return np.array(data_diff)
  

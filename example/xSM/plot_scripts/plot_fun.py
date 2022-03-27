import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def fun_gamma(data):
  np.seterr(divide='ignore',invalid='ignore')
  return abs(data[:,5]-data[:,7])/data[:,4]

def fun_gamma_line(line):
  return abs(line[5]-line[7])/line[4]
  
def loaddata(name, gamma_min=0.1):
  data = np.loadtxt(name)
  data = data[data[:,3]>0]
  return data[fun_gamma(data)>gamma_min]

def fun_diff(data1, data2, data0, show_gamma=False, norm =False, use_abs=False, gamma_min = 0, sort=True):
  if len(data1) != len(data2) or len(data1) != len(data0):
    print("Wrong data files")
    sys.exit()
  data_diff = []
  for ii in range(len(data1)):
    if data1[ii][0] != data2[ii][0] or data1[ii][0] != data0[ii][0]:
      print("Wrong data files")
      print(ii, data1[ii][0], data2[ii][0], data0[ii][0])
      sys.exit()
    ms = data1[ii][0]
    lambda_s = data1[ii][1]
    lambda_hs = data1[ii][2]

    size_1 = data1[ii][3]
    size_2 = data2[ii][3]
    size_0 = data0[ii][3]
    
    if size_1 > 0 and size_2 > 0 and size_0>0:
      TC_1 = data1[ii][4]
      TC_2 = data2[ii][4]
      TC_0 = data0[ii][4]
      vs_0 = data0[ii][8]
      gamma_1 = fun_gamma_line(data1[ii])
      gamma_2 = fun_gamma_line(data2[ii])
      gamma_0 = fun_gamma_line(data0[ii])
      if gamma_0 >= gamma_min and data1[ii][8]>10 and data2[ii][8]>10:
        if show_gamma:
          diff = abs(gamma_1 - gamma_2) if use_abs else gamma_1 - gamma_2
          x1 = gamma_1
          x2 = gamma_2
          central = gamma_0
        else:
          diff = abs(TC_1 - TC_2) if use_abs else TC_1 - TC_2
          x1 = TC_1
          x2 = TC_2
          central = TC_0
#          if diff > 20:
#            print "---------------------"
#            print "input = ", data1[ii][0:3]
#            print "T = ", x1 ,x2
#            print "vh = ", data1[ii][5], data2[ii][5]
#            print "vs = ", data1[ii][8], data2[ii][8]
            
             
        data_diff.append([ms, lambda_s, lambda_hs, 1, diff/central if norm else diff, central, x1, x2])
        
  if sort:
    data_diff.sort(key=(lambda x:abs(x[4]))) 
  
  return np.array(data_diff)


#def fun_diff(data1, data2, data0, show_gamma=False, norm =False, use_abs=False, gamma_min = 0.4, sort=True):
#  if len(data1) != len(data2) or len(data1) != len(data0):
#    print("Wrong data files")
#    sys.exit()
#  data_diff = []
#  data_1 = []
#  data_2 = []
#  for ii in range(len(data1)):
#    if data1[ii][0] != data2[ii][0] or data1[ii][0] != data0[ii][0]:
#      print("Wrong data files")
#      print(ii, data1[ii][0], data2[ii][0], data0[ii][0])
#      sys.exit()
#    ms = data1[ii][0]
#    lambda_s = data1[ii][1]
#    lambda_hs = data1[ii][2]

#    size_1 = data1[ii][3]
#    size_2 = data2[ii][3]
#    size_0 = data0[ii][3]
#    
#    if size_1 > 0 and size_2 > 0 and size_0>0:
#      TC_1 = data1[ii][4]
#      TC_2 = data2[ii][4]
#      TC_0 = data0[ii][4]
#      vs_0 = data0[ii][8]
#      if TC_1>0 and TC_0>0 and TC_0>0  \
#        and abs(data0[ii][5])>1 and abs(data0[ii][6])<1 and abs(data0[ii][7])<1 and abs(data0[ii][8])>1:
#        gamma_1 = fun_gamma_line(data1[ii])
#        gamma_2 = fun_gamma_line(data2[ii])
#        gamma_0 = fun_gamma_line(data0[ii])
#        if gamma_0 >= gamma_min and gamma_0 < 10: 
#          if show_gamma:
#            diff = abs(gamma_1 - gamma_2) if use_abs else gamma_1 - gamma_2
#            x1 = gamma_1
#            x2 = gamma_2
#            central = gamma_0
#          else:
#            diff = abs(TC_1 - TC_2) if use_abs else TC_1 - TC_2
#            x1 = TC_1
#            x2 = TC_2
#            central = TC_0
#          data_diff.append([ms, lambda_s, lambda_hs, diff/central if norm else diff, central, x1, x2])
#    elif size_1 != 0 and size_2 == 0:
#      if gamma_1 > gamma_min :
#        data_1.append([lambda_hs, lambda_s, ms])
#    elif size_1 == 0 and size_2 != 0:
#      if gamma_2 > gamma_min :
#        data_2.append([lambda_hs, lambda_s, ms])
#        
#        
#  if sort:
#    data_diff.sort(key=(lambda x:abs(x[3]))) 
#  
#  return np.array(data_diff), np.array(data_1), np.array(data_2)   
  

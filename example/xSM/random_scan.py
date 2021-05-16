import os,sys,signal
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()
print cwd

n_total = 100000
#n_total = 100

folder_name = "random_scan_results"
if not os.path.exists(folder_name):
  os.mkdir(folder_name) 

def scan(cmd, file_name):
  os.chdir(cwd)
  fo = open(folder_name+"/"+file_name+".txt", "w")
  if not os.path.exists(file_name):
    os.mkdir(file_name)
  os.chdir(cwd+"/"+file_name)
  for ii in range(n_total):
    par = (" " + str(ms[ii]) + " " +str(lambda_s[ii]) + " " + str(lambda_hs[ii]) 
          +" " + str(Q) + " " + str(xi) + " " + str(daisy_flag) 
          +" " + str(use_1L_EWSB_in_0L_mass) )

    print par
    os.system(cmd+par)
    output = open("output.txt").readline()
    fo.write( output )
  fo.close()


# for run_ScalarSingletZ2DMMhInput_withSingletVEVinPT, ms means mu_s
np.random.seed(1)
ms = np.random.uniform(10,110,n_total)
lambda_s = np.random.uniform(0.01,0.3,n_total)
lambda_hs = np.random.uniform(0.1,0.5,n_total)

Q = 173
xi = 1
daisy_flag = 1
use_1L_EWSB_in_0L_mass = 0

######################################
#scheme = "xSM_MSbar"
#cmd = "./../../../bin/run_"+scheme

#xi = 1
#file_name = scheme + "_xi_1"
#scan(cmd, file_name)

#xi = 0.1
#file_name = scheme + "_xi_01"
#scan(cmd, file_name)

#xi = 3
#file_name = scheme + "_xi_3"
#scan(cmd, file_name)

######################################
#daisy_flag = 1
#file_name = scheme + "_daisy_Parwani"
#scan(cmd, file_name)

#daisy_flag = 2
#file_name = scheme + "_daisy_ArnoldEspinosa"
#scan(cmd, file_name)


#####################################
scheme = "ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
cmd = "./../../../bin/run_"+scheme

#file_name = scheme + "_mt"
#Q = 173.
#scan(cmd, file_name)

#file_name = scheme + "_2mt"
#Q = 173.*2.
#scan(cmd, file_name)

file_name = scheme + "_05mt"
Q = 173./2.
scan(cmd, file_name)



import os,sys,signal
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()
print cwd

n_total = 100000
#n_total = 10

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
          +" " + str(use_1L_EWSB_in_0L_mass) + " " + str(use_Goldstone_resum) )

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

mt = 173.
Q = mt
xi = 1
daisy_flag = 2
use_1L_EWSB_in_0L_mass = 0
use_Goldstone_resum = 0


####################################
#scheme = "ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#cmd = "./../../../bin/run_"+scheme

#file_name = scheme + "_mt"
#Q = mt
#scan(cmd, file_name)

#file_name = scheme + "_2mt"
#Q = mt*2.
#scan(cmd, file_name)

#file_name = scheme + "_05mt"
#Q = 173./2.
#scan(cmd, file_name)


#####################################
scheme = "xSM_MSbar"
cmd = "./../../../bin/run_xSM_MSbar"
use_1L_EWSB_in_0L_mass = 0
use_Goldstone_resum = 1

#file_name = scheme+"xi0"
#xi = 0
#scan(cmd, file_name)

file_name = scheme+"xi25"
xi = 25
scan(cmd, file_name)




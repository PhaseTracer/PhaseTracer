import os,sys,signal
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()
print cwd

folder_name = "scan_results"
if not os.path.exists(folder_name):
  os.mkdir(folder_name) 

def scan(cmd, file_name):
  os.chdir(cwd)
  fo = open(folder_name+"/"+file_name+".txt", "w")
  if not os.path.exists(file_name):
    os.mkdir(file_name)
  os.chdir(cwd+"/"+file_name)
  for ms_ in ms:
    for lambda_s_ in lambda_s:
      for lambda_hs_ in lambda_hs:
        for Q_ in Q:
          for xi_ in xi:
            for daisy_flag_ in daisy_flag:
              for use_1L_EWSB_in_0L_mass_ in use_1L_EWSB_in_0L_mass:

                par = (" " + str(ms_) + " " +str(lambda_s_) + " " + str(lambda_hs_) 
                      +" " + str(Q_) + " " + str(xi_) + " " + str(daisy_flag_) 
                      +" " + str(use_1L_EWSB_in_0L_mass_) )

                print par
                os.system(cmd+par)
                os.system(cmd+par)
                output = open("output.txt").readline()
                fo.write( output )
  fo.close()


scan_ms = True
scan_lambda_s = True
scan_lambda_hs = True

# for run_ScalarSingletZ2DMMhInput_withSingletVEVinPT, ms means mu_s
ms = np.linspace(10,110,50) if scan_ms else [62.5]
lambda_s = np.linspace(0.01,0.3,50) if scan_lambda_s else [0.1]
lambda_hs = np.linspace(0.1,0.5,50) if scan_lambda_hs else [0.25]

Q = [173]
xi = [1]
daisy_flag = [1]
use_1L_EWSB_in_0L_mass = [0]


scheme = "xSM_MSbar"
cmd = "./../../../bin/run_"+scheme

#xi = [1]
#file_name = scheme + "_xi_1"
#scan(cmd, file_name)

#xi = [0.1]
#file_name = scheme + "_xi_01"
#scan(cmd, file_name)

#xi = [3]
#file_name = scheme + "_xi_3"
#scan(cmd, file_name)

######################################
#daisy_flag = [1]
#file_name = scheme + "_daisy_Parwani"
#scan(cmd, file_name)

daisy_flag = [2]
file_name = scheme + "_daisy_ArnoldEspinosa"
scan(cmd, file_name)


######################################
#scheme = "xSM_MSbar"
#cmd = "./../../../bin/run_"+scheme

#file_name = scheme + "_mt"
#Q = [173.]
#scan(cmd, file_name)

#file_name = scheme + "_2mt"
#Q = [173.*2.]
#scan(cmd, file_name)

#file_name = scheme + "_05mt"
#Q = [173./2.]
#scan(cmd, file_name)
#Q = [173.]


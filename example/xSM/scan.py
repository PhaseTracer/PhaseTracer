import os,sys,signal
import re,shutil
import subprocess
import numpy as np
import time


scan_ms = True
scan_lambda_s = False
scan_lambda_hs = False

# for run_ScalarSingletZ2DMMhInput_withSingletVEVinPT, ms means mu_s
ms = np.linspace(30,60,5) if scan_ms else [62.5]
lambda_s = np.linspace(0.01,0.2,5) if scan_lambda_s else [0.1]
lambda_hs = np.linspace(0.1,0.4,5) if scan_lambda_hs else [0.25]
Q = [173]
xi = [1]
daisy_flag = [1]
use_1L_EWSB_in_0L_mass = [0]


file_name = "scan_results"
if not os.path.exists(file_name):
  os.mkdir(file_name) 
file_name += "/"

def scan(cmd, file_name):
  if scan_ms: file_name += "_ms"
  if scan_lambda_s: file_name += "_lambda_s"
  if scan_lambda_hs: file_name += "_lambda_hs"
  
  fo = open(file_name+".txt", "w")
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
                
                output = open("output.txt").readline()
                
                fo.write( output )
  fo.close()


scheme = "xSM_MSbar"
file_name += scheme
cmd = "./../../bin/run_"+scheme
scan(cmd, file_name)







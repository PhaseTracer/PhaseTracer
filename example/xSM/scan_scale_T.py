import os,sys,signal,shutil
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()

scan_1d_bks = True
scan_ms = False
scan_ls = False
scan_lhs = True




if scan_1d_bks or scan_1d_xi:
  n_total = 200
  folder_name = "1d_bks"


if not os.path.exists(folder_name):
  os.mkdir(folder_name) 

pi2 = 3.1415926*2

def runPT(Q_, other_settings):
  cmd, ms_, lambda_s_, lambda_hs_, xi_, daisy_flag, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, special_flag = other_settings
  par = (" " + str(ms_) + " " +str(lambda_s_) + " " + str(lambda_hs_) 
        +" " + str(Q_) + " " + str(xi_) + " " + daisy_flag 
        +" " + use_1L_EWSB_in_0L_mass + " "+ use_Goldstone_resum 
        +" " + special_flag)
  os.system(cmd+par)
  
def iterate(Q_, k, other_settings):
  Q_in = Q_
  for i in range(100):
    runPT(k*Q_in, other_settings)
    outpar = np.loadtxt("output.txt")
    if len(outpar)< 5:
      return  Q_in/k
    Q_out = outpar[4]
    print(Q_in,Q_out)
    if abs(Q_out - Q_in) < 0.1:
      return Q_out
    Q_in = Q_out
  return Q_out
                
def scan(cmd, file_name, ms, lambda_s, lambda_hs, xi, daisy_flag, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, special_flag=""):
  os.chdir(cwd)
  fo_T = open(folder_name+"/"+file_name+"T.txt", "w")
  fo_2piT = open(folder_name+"/"+file_name+"2piT.txt", "w")
  fo_05T = open(folder_name+"/"+file_name+"05T.txt", "w")
  if not os.path.exists(file_name):
    os.mkdir(file_name)
  os.chdir(cwd+"/"+file_name)
  for ms_ in ms:
    for lambda_s_ in lambda_s:
      for lambda_hs_ in lambda_hs:
          for xi_ in xi:
                other_settings = [cmd, ms_, lambda_s_, lambda_hs_, xi_, daisy_flag, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, special_flag]
                runPT(173, other_settings)
                T_mt = np.loadtxt("output.txt")[4]
                if (T_mt>5 and T_mt < 200 ):

                  T_out = iterate(T_mt, 1., other_settings)
                  runPT(T_out, other_settings)
                  output = open("output.txt").readline()
                  fo_T.write( output )

                  T_out = iterate(T_mt, pi2, other_settings)
                  runPT(pi2*T_out, other_settings)
                  output = open("output.txt").readline()
                  fo_2piT.write( output )
                  
                  T_out = iterate(T_mt, 0.5, other_settings)
                  runPT(0.5*T_out, other_settings)
                  output = open("output.txt").readline()
                  fo_05T.write( output )
                  
#                  runPT(T_mt, other_settings)
#                  output = open("output.txt").readline()
#                  fo_T.write( output )

#                  runPT(pi2*T_mt, other_settings)
#                  output = open("output.txt").readline()
#                  fo_2piT.write( output )
#                  
#                  runPT(0.5*T_mt, other_settings)
#                  output = open("output.txt").readline()
#                  fo_05T.write( output )
                  
                  
                  
  fo_T.close()
  fo_2piT.close()
  fo_05T.close()
  os.chdir(cwd)
  shutil.rmtree(file_name)

# default choice
xi_in = [1]
daisy_flag_in = "2" # ArnoldEspinosa
use_1L_EWSB_in_0L_mass_in = "0"
use_Goldstone_resum_in = "1"

def perfrom_1d_scan(ms, lambda_s, lambda_hs, file_name_):

#  ############### scale band ####################

  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
  file_name = file_name_
  scan(cmd, file_name, ms, lambda_s, lambda_hs, [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"PW_"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [1], "1", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noRGE_woFS_"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "1")

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_noRGE_"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in,"1")



if scan_1d_bks:
  if scan_lhs:
    #1d scan for lambda_hs
    ms = [65]
    lambda_s = [0.1]
    lambda_hs = np.linspace(0.1,0.4,n_total)
    filename = "lambda_hs_"
      
    perfrom_1d_scan(ms, lambda_s, lambda_hs, filename)
    
  if scan_ls:
    #1d scan for lambda_s
    ms = [65]
    lambda_hs = [0.3]
    lambda_s = np.linspace(0.01,0.2,n_total)
    filename = "lambda_s_"
      
    perfrom_1d_scan(ms, lambda_s, lambda_hs, filename)
  
  if scan_ms:
    #1d scan for m_s
    lambda_s = [0.1]
    lambda_hs = [0.3]
    ms = np.linspace(40,110,n_total)
    filename = "m_s_"
    perfrom_1d_scan(ms, lambda_s, lambda_hs, filename)


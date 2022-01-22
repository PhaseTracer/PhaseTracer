import os,sys,signal,shutil
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()

scan_1d_bks = True
scan_ms = True
scan_ls = True
scan_lhs = True
scan_1d_xi = False
xi_zoom_in = False


scan_2d_scan = False 
scan_ls_lhs = False
scan_ms_ls = False
scan_ms_lhs = False



if scan_1d_bks or scan_1d_xi:
  n_total = 200
  folder_name = "1d_bks"

if scan_2d_scan:
  n_x = 200
  n_y = n_x
  folder_name = "2d_scan"



if not os.path.exists(folder_name):
  os.mkdir(folder_name) 

def scan(cmd, file_name, ms, lambda_s, lambda_hs, Q, xi, daisy_flag, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, special_flag=""):
  print(xi)
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
                par = (" " + str(ms_) + " " +str(lambda_s_) + " " + str(lambda_hs_) 
                      +" " + str(Q_) + " " + str(xi_) + " " + daisy_flag 
                      +" " + use_1L_EWSB_in_0L_mass + " "+ use_Goldstone_resum 
                      +" " + special_flag)

                print(par)
                os.system(cmd+par)
                output = open("output.txt").readline()
                fo.write( output )
  fo.close()
  os.chdir(cwd)
  shutil.rmtree(file_name)

# default choice
Q_in = [173]
xi_in = [0]
daisy_flag_in = "2" # ArnoldEspinosa
use_1L_EWSB_in_0L_mass_in = "0"
use_Goldstone_resum_in = "1"

def perfrom_1d_scan(ms, lambda_s, lambda_hs, file_name_):

###  ############## methods ################
#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"default"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_0L"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "1")

#  cmd = "./../../../bin/run_xSM_OSlike"
#  file_name = file_name_+"OSlike"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"Parwani"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, "1", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noDaisy"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_HT"
#  file_name = file_name_+"HT"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)


#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"covariant_gauge"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "2")

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"xi3"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [3], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"1L_EWSB"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, "1", "0")

#  ############### scale band ####################

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "1")
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "1")
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "1")


#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_noRGE_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "3")
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_noRGE_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "3")
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_0L_noRGE_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "3")

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173*2], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noRGE_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "3")

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noRGE_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173*2], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "3")

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noRGE_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "3")


#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173*2], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noRGE_woFS_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noRGE_woFS_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noRGE_woFS_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)  

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_xi0"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_noRGE_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in,"1")

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_noRGE_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in,"1")
#  
#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_noRGE_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in,"1")

#  ############### xi band ####################

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"xi3"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [3], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"xi1"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"xi0"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noD_xi3"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [3], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noD_xi1"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"noD_xi0"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_xi3"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [3], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_xi1"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [1], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)
#  
#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_xi0"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

  cmd = "./../../../bin/run_xSM_PRM"
  file_name = file_name_+"PRM_woFS_0L"
  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], "0", use_1L_EWSB_in_0L_mass_in, "0", "1")

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
    lambda_s = np.linspace(0.01,0.2,n_total)
    lambda_hs = [0.3]
    filename = "lambda_s_"
    perfrom_1d_scan(ms, lambda_s, lambda_hs, filename)
  
  if scan_ms:
    #1d scan for m_s
    ms = np.linspace(40,110,n_total)
    lambda_s = [0.1]
    lambda_hs = [0.3]
    filename = "m_s_"
    perfrom_1d_scan(ms, lambda_s, lambda_hs, filename)
  
  if scan_1d_xi:
    ms = [65]
    lambda_s = [0.1]
    lambda_hs = [0.3]
    if xi_zoom_in:
      xi = np.linspace(0,0.2,n_total)
      add_file_name = "_zoom_in"
    else:
      xi = np.linspace(0,10,n_total)
      add_file_name = ""
      
    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "Rxi_MSbar"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "covariant_MSbar"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "2")

    cmd = "./../../../bin/run_xSM_HT"
    file_name = "Rxi_HT"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "Rxi_PRM"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "Rxi_PRM_0L"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, "0", "1")

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "Rxi_PRM_0L_resum"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, "1", "1")

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "covariant_PRM"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", "0", "1", "2")

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "covariant_PRM_0L"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", "0", "0", "3")

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "covariant_PRM_0L_resum"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", "0", "1", "3")

    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "Rxi_MSbar"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

##    xi = np.linspace(0,2,n_total)
##    cmd = "./../../../bin/run_xSM_MSbar"
##    file_name = "Rxi_MSbar_resummation"
##    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)
    
    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "Rxi_MSbar_1L_EWSB"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, "1", "0")
 
    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "Rxi_MSbar_no"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, "0", "0") 

    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "covariant_MSbar_1L_EWSB"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, "1", "0", "2")

    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "covariant_MSbar_no"+add_file_name
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, "0", "0", "2")










def perfrom_2d_scan(ms, lambda_s, lambda_hs, file_name_):

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"default"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"xi3"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [3], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
  file_name = file_name_+"05mt"
  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], [1], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_OSlike"
#  file_name = file_name_+"OSlike"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)


#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"Parwani"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, "1", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)


if scan_2d_scan:
  
  if scan_ls_lhs:
    ms = [65]
    lambda_s = np.linspace(0.01,0.3,n_x)
    lambda_hs = np.linspace(0.1,0.5,n_y)
    filename = "lhs_ls_"
    perfrom_2d_scan(ms, lambda_s, lambda_hs, filename)

  if scan_ms_ls:
    ms = np.linspace(10,120,n_x)
    lambda_s = np.linspace(0.01,0.3,n_y)
    lambda_hs = [0.3]
    filename = "ms_ls_"
    perfrom_2d_scan(ms, lambda_s, lambda_hs, filename)
    
  if scan_ms_lhs:
    ms = np.linspace(10,120,n_x)
    lambda_s = [0.1]
    lambda_hs = np.linspace(0.1,0.5,n_y)
    filename = "ms_lhs_"
    perfrom_2d_scan(ms, lambda_s, lambda_hs, filename)



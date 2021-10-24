import os,sys,signal,shutil
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()
print cwd


scan_1d_bks = True
scan_ms = True
scan_ls = True
scan_lhs = True
scan_1d_xi = False

scan_2d_bks = True


if scan_1d_bks or scan_1d_xi:
  n_total = 100
  folder_name = "1d_bks"


if not os.path.exists(folder_name):
  os.mkdir(folder_name) 

def scan(cmd, file_name, ms, lambda_s, lambda_hs, Q, xi, daisy_flag, use_1L_EWSB_in_0L_mass, use_Goldstone_resum, special_flag=""):
  print xi
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

                print par
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

#  ############## methods ################
#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"default"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_xSM_OSlike"
#  file_name = file_name_+"OSlike"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)


#  cmd = "./../../../bin/run_xSM_HT"
#  file_name = file_name_+"HT"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)


#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"covariant_gauge"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "2")

#  cmd = "./../../../bin/run_xSM_MSbar"
#  file_name = file_name_+"xi3"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [3], daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

  cmd = "./../../../bin/run_xSM_MSbar"
  file_name = file_name_+"1L_EWSB"
  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi_in, daisy_flag_in, "1", "0")

#  ############### scale band ####################

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)
#  
#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT_PRM"
#  file_name = file_name_+"PRM_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [2*173], xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173*2], xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], xi_in, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173], xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_2mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [173*2], xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

#  cmd = "./../../../bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#  file_name = file_name_+"noD_05mt"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, [0.5*173], xi_in, "0", use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

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

#  cmd = "./../../../bin/run_xSM_PRM"
#  file_name = file_name_+"PRM_woFS_0L"
#  scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, [0], "0", use_1L_EWSB_in_0L_mass_in, "0", "1")

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
    xi = np.linspace(0,10,n_total)

    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "Rxi_MSbar"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

    cmd = "./../../../bin/run_xSM_MSbar"
    file_name = "covariant_MSbar"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in, "2")

    cmd = "./../../../bin/run_xSM_HT"
    file_name = "Rxi_HT"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, daisy_flag_in, use_1L_EWSB_in_0L_mass_in, use_Goldstone_resum_in)

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "Rxi_PRM"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, "0")

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "Rxi_PRM_0L"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, "0", "1")
    
    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "covariant_PRM"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, "0", "2")

    cmd = "./../../../bin/run_xSM_PRM"
    file_name = "covariant_PRM_0L"
    scan(cmd, file_name, ms, lambda_s, lambda_hs, Q_in, xi, "0", use_1L_EWSB_in_0L_mass_in, "0", "3")
    
    
    
#scheme = "xSM_MSbar"
#cmd = "./../../../bin/run_"+scheme
##################### gauge dependence #################
# 2d 
#scheme = "xSM_MSbar"
#cmd = "./../../../bin/run_"+scheme
#n_total = 200
#
#xi = [0]
#file_name_ = scheme + "_xi_0"
#
#file_name  = file_name_+ "_ms_fixed"
#ms = [65]
#lambda_s = np.linspace(0.01,0.3,n_total)
#lambda_hs = np.linspace(0.1,0.5,n_total)
#scan(cmd, file_name)
#
#
#file_name  = file_name_+ "_ls_fixed"
#ms = np.linspace(10,110,n_total)
#lambda_s = [0.1]
#lambda_hs = np.linspace(0.1,0.5,n_total)
#scan(cmd, file_name)
#
#
#file_name  = file_name_+ "_lhs_fixed"
#ms = np.linspace(10,110,n_total)
#lambda_s = np.linspace(0.01,0.3,n_total)
#lambda_hs = [0.3]
#scan(cmd, file_name)


##################### gauge dependence #################
#BK = 'BK3_0'
#ms =  [75.4206]
#lambda_s =  [0.0390329]
#lambda_hs =  [0.32759]
#
#BK = 'BK3_ms_p'
#ms =  [77]
#lambda_s =  [0.0390329]
#lambda_hs =  [0.32759]
#
#BK = 'BK3_ms_m'
#ms =  [75]
#lambda_s =  [0.0390329]
#lambda_hs =  [0.32759]
##
#BK = 'BK3_ls_p'
#ms =  [75.4206]
#lambda_s =  [0.041]
#lambda_hs =  [0.32759]
#
#BK = 'BK3_ls_m'
#ms =  [75.4206]
#lambda_s =  [0.037]
#lambda_hs =  [0.32759]
#
#BK = 'BK3_lhs_p'
#ms =  [75.4206]
#lambda_s =  [0.0390329]
#lambda_hs =  [0.33]

#BK = 'BK3_lhs_m'
#ms =  [75.4206]
#lambda_s =  [0.0390329]
#lambda_hs =  [0.29]

#BK = 'BK1'
#ms =  [106.918]
#lambda_s =  [0.0190792]
#lambda_hs =  [0.474493]
#BK = 'BK2'
#ms =  [85.4066]
#lambda_s =  [0.0130036]
#lambda_hs =  [0.319915]
##BK = 'BK3'
##ms =  [75.4206]
##lambda_s =  [0.0390329]
##lambda_hs =  [0.32759]
##BK = 'BK4'
##ms =  [105.836]
##lambda_s =  [0.0304002]
##lambda_hs =  [0.446975]
##BK = 'BK5'
##ms =  [103.036]
##lambda_s =  [0.0654765]
##lambda_hs =  [0.41981]
##BK = 'BK6'
##ms =  [18.5591]
##lambda_s =  [0.271993]
##lambda_hs =  [0.387765]

#file_name = BK

#xi = np.linspace(0, 5, 30)

#file_name = BK + "_goldstone_resum"
#use_1L_EWSB_in_0L_mass = "0"
#use_Goldstone_resum = "1"


##file_name = BK + "_use_1L_EWSB_in_0L_mass"
##use_1L_EWSB_in_0L_mass = "1"
##use_Goldstone_resum = "0"

##file_name = BK + "_no_gildstone_catastrophe"
##use_1L_EWSB_in_0L_mass = "0"
##use_Goldstone_resum = "0"

#file_name = BK + "_both"
#use_1L_EWSB_in_0L_mass = "1"
#use_Goldstone_resum = "1"

#scan(cmd, file_name)

######################################
#daisy_flag = "1"
#file_name = scheme + "_daisy_Parwani"
#scan(cmd, file_name)

#daisy_flag = "2"
#file_name = scheme + "_daisy_ArnoldEspinosa"
#scan(cmd, file_name)


######################################
#scheme = "ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
#cmd = "./../../../bin/run_"+scheme
#mt = 173.
#use_Goldstone_resum = "0"

#BK = 'BK1'
#ms =  [37.1949]
#lambda_s =  [0.281362]
#lambda_hs =  [0.422943]
#BK = 'BK2'
#ms =  [70.4449]
#lambda_s =  [0.013275]
#lambda_hs =  [0.199424]
#BK = 'BK3'
#ms =  [82.6337]
#lambda_s =  [0.0595215]
#lambda_hs =  [0.291194]
#BK = 'BK4'
#ms =  [19.9661]
#lambda_s =  [0.149159]
#lambda_hs =  [0.290839]
#BK = 'BK5'
#ms =  [93.4236]
#lambda_s =  [0.0206235]
#lambda_hs =  [0.457805]
#BK = 'BK6'
#ms =  [42.3438]
#lambda_s =  [0.113348]
#lambda_hs =  [0.286329]

#file_name = BK + "_mu"
#Q = np.linspace(0.5*mt, 2*mt, 20)
#scan(cmd, file_name)



#file_name = scheme + "_2mt"
#Q = [173.*2.]
#scan(cmd, file_name)

#file_name = scheme + "_05mt"
#Q = [173./2.]
#scan(cmd, file_name)
#Q = [173.]


###################### gauge dependence #################

#scheme = "xSM_OSlike"
#cmd = "./../../../bin/run_"+scheme
#n_total = 50

#xi = [0]
#file_name_ = scheme + "_xi_0"

#file_name  = file_name_+ "_ms_fixed"
#ms = [65]
#lambda_s = np.linspace(0.01,0.3,n_total)
#lambda_hs = np.linspace(0.1,0.5,n_total)
#scan(cmd, file_name)


#file_name  = file_name_+ "_ls_fixed"
#ms = np.linspace(10,110,n_total)
#lambda_s = [0.1]
#lambda_hs = np.linspace(0.1,0.5,n_total)
#scan(cmd, file_name)


#file_name  = file_name_+ "_lhs_fixed"
#ms = np.linspace(10,110,n_total)
#lambda_s = np.linspace(0.01,0.3,n_total)
#lambda_hs = [0.3]
#scan(cmd, file_name)


#scheme = "xSM_OSlike"
#cmd = "./../../../bin/run_"+scheme

#file_name = "BK1"
#xi = np.linspace(0, 5, 30)
#ms =  [37.1949]
#lambda_s =  [0.281362]
#lambda_hs =  [0.422943]
#
#file_name = "BK1_no_gildstone_catastrophe"
#xi = np.linspace(0.01, 5, 30)
#use_Goldstone_resum = "0"
#ms =  [37.1949]
#lambda_s =  [0.281362]
#lambda_hs =  [0.422943]
#
#scan(cmd, file_name)

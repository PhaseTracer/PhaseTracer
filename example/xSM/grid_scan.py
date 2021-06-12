import os,sys,signal
import re,shutil
import subprocess
import numpy as np
import time

cwd = os.getcwd()
print cwd


#folder_name = "scan_results"

#folder_name = "gauge_dependence"

folder_name = "goldstone_catastrophe"

#folder_name = "scale_dependence"

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
                par = (" " + str(ms_) + " " +str(lambda_s_) + " " + str(lambda_hs_) 
                      +" " + str(Q_) + " " + str(xi_) + " " + daisy_flag 
                      +" " + use_1L_EWSB_in_0L_mass + " "+ use_Goldstone_resum )

                print par
                os.system(cmd+par)
                output = open("output.txt").readline()
                fo.write( output )
  fo.close()


scan_ms = True
scan_lambda_s = True
scan_lambda_hs = True

Q = [173]
xi = [1]
daisy_flag = "2" # ArnoldEspinosa
use_1L_EWSB_in_0L_mass = "0"
use_Goldstone_resum = "1"

#scheme = "xSM_MSbar"
#cmd = "./../../../bin/run_"+scheme
##################### gauge dependence #################

##BK = 'BK1'
##ms =  [106.918]
##lambda_s =  [0.0190792]
##lambda_hs =  [0.474493]
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

#xi = np.linspace(0, 5, 200)

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
scheme = "ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT"
cmd = "./../../../bin/run_"+scheme
mt = 173.
use_Goldstone_resum = "0"

#BK = 'BK1'
#ms =  [37.1949]
#lambda_s =  [0.281362]
#lambda_hs =  [0.422943]
#BK = 'BK2'
#ms =  [70.4449]
#lambda_s =  [0.013275]
#lambda_hs =  [0.199424]
BK = 'BK3'
ms =  [82.6337]
lambda_s =  [0.0595215]
lambda_hs =  [0.291194]
BK = 'BK4'
ms =  [19.9661]
lambda_s =  [0.149159]
lambda_hs =  [0.290839]
BK = 'BK5'
ms =  [93.4236]
lambda_s =  [0.0206235]
lambda_hs =  [0.457805]
#BK = 'BK6'
#ms =  [42.3438]
#lambda_s =  [0.113348]
#lambda_hs =  [0.286329]

file_name = BK + "_mu"
Q = np.linspace(0.5*mt, 2*mt, 20)
scan(cmd, file_name)



#file_name = scheme + "_2mt"
#Q = [173.*2.]
#scan(cmd, file_name)

#file_name = scheme + "_05mt"
#Q = [173./2.]
#scan(cmd, file_name)
#Q = [173.]



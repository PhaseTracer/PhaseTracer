### Scan

    python grid_scan.py


  The scan result are saved in scan_results folder. 

  xSM_MSbar_xi_01.txt means xi=0.1, using Parwani method, Q = 173
  xSM_MSbar_xi_1.txt means xi=1, other setting are same as above
  xSM_MSbar_xi_3.txt means xi=3, other setting are same as above
  
  xSM_MSbar_daisy_ArnoldEspinosa.txt means using ArnoldEspinosa method, xi=1, Q = 173
  xSM_MSbar_daisy_Parwani.txt means using Parwani method,  other setting are same as above
  
  
  The meanings of each column are
    0 : ms 
    1 : lambda_s
    2 : lambda_hs
    3 : Number of transitions. =-1 means can not find phase; =-2 means can not find any transition
    4 : TC of the transition with largest gamma
    5 : true_vacuum[0] of the transition with largest gamma
    6 : true_vacuum[1] of the transition with largest gamma
    7 : false_vacuum[0] of the transition with largest gamma
    8 : false_vacuum[1] of the transition with largest gamma
    9 : Renormalization scale
    10 : xi, gauge fixing parameter
    11 : daisy_flag. =0 means no daisy term; =1 means Parwani; =2 means ArnoldEspinosa
    12: use_1L_EWSB_in_0L_mass. =0 means use resum Goldstone contribution; =1 means use_1L_EWSB_in_0L_mass

### Plot

    python compare.py
    
  There are some flags you can change in compare.py.

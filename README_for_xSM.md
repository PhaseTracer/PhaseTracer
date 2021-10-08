Basically, I use one cpp main file for each scheme, and the c++ codes only calculate one point in parameter space. I use python scripts to scan the parameter space by calling the c++ codes, so that it is easy to change the scan method and ranges. 

### Building

    mkdir build
    cd build
    cmake ..
    make
    
 To build the one using FS to run RGE, 
    cmake -DBUILD_WITH_FS_ScalarSingletZ2DMMhInput=on ..
    make
    
### run one point

  ./bin/run_xSM_MSbar x1 x2 x3 x4 x5 x6 x7
  ./bin/run_ScalarSingletZ2DMMhInputMsInput_withSingletVEVinPT x1 x2 x3 x4 x5 x6 x7
  ./bin/run_xSM_OSlike x1 x2 x3 x4 x5 x6 x7
  ./bin/run_xSM_PRM x1 x2 x3 x4 x5 x6 x7
  ./bin/run_xSM_HT x1 x2 x3

  where
  x1: ms
  x2: lambda_s
  x3: lambda_hs
  x4: renormalization scale, without RGE running
  x5: xi, gauge fixing parameter
  x6: daisy_flag: 0 -> no daisy correction
                  1 -> Parwani
                  2 -> ArnoldEspinosa
  x7: use_1L_EWSB_in_0L_mass: 0 -> false
                              1 -> true
                       
  
  such as 
  ./bin/run_xSM_MSbar 60 0.1 0.25 173 1 1 0
                       
  Outputs are in output.txt, in order of 
    ms, lambda_s, lambda_hs, number_of_transitions, TC, true_vacuum[0], true_vacuum[1], false_vacuum[0], false_vacuum[1], Q, xi, daisy_flag, use_1L_EWSB_in_0L_mass
  
  where number_of_number_of_transitionss = -1 or -2 means that can not find phase or transition for this point. 
  I only saved Tc, vh, vs of the strongest transition. 
 
 
### scan

  cd example/xSM/
  python grid_scan.py

  There are many flags in grid_scan.py. Turn on them to perfrom scans you need.


  In the output files, the meanings of each column are same to the above output file, which are
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


### plot

  cd example/xSM/plot_scripts/
  python 1d_plot_bks.py
  python 1d_plot_scale.py
  







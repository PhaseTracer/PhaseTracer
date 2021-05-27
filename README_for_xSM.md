### Building

    mkdir build
    cd build
    cmake ..
    make
    
 To build the one using FS to run RGE, 
    cmake -DBUILD_WITH_FS_ScalarSingletZ2DMMhInput=on ..
    make
    
### run MSbar scheme without RGE running
  ./bin/run_xSM_MSbar x1 x2 x3 x4 x5 x6 x7
  
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
                              
### run MSbar scheme with RGE running
  /bin/run_ScalarSingletZ2DMMhInput_withSingletVEVinPT x1 x2 x3 x4 x5 x6 x7
  
  where
  x1: mu_s^2
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
  ./bin/run_ScalarSingletZ2DMMhInput_withSingletVEVinPT -3000 0.1 0.25 173 1 1 0    
                      
  Output is same to the run_xSM_MSbar except that "ms" is mu_s^2
  
### run Onshell like scheme
        
  /bin/run_xSM_OSlike x1 x2 x3 x4 x5 x6
  
  where
  x1: mu_s^2
  x2: lambda_s
  x3: lambda_hs
  x4: xi, gauge fixing parameter
  x5: daisy_flag: 0 -> no daisy correction
                  1 -> Parwani
                  2 -> ArnoldEspinosa
  x6: use_Goldstone_resum: 0 -> false  # Must be true if xi=0
                           1 -> true
                           
  such as                         
  .bin/run_xSM_OSlike 62.5 0.124763 0.24 0 1 1
  
  Outputs are in output.txt, in order of 
    ms, lambda_s, lambda_hs, number_of_transitions, TC, true_vacuum[0], true_vacuum[1], false_vacuum[0], false_vacuum[1], xi, daisy_flag, use_Goldstone_resum
    
    
                  
                              

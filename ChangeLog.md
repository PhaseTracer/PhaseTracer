## PhaseTracer-1.0.0 [March 4, 2020]
 * Initial release
## PhaseTracer-1.0.1 [April 10, 2020]
 * Setup automatic build & tests etc
 * BugFix: calculate_sm_masses setting in FlexibleSUSY should be set to 1 or 0.
 * Fix cmake building of FS example to always use tagged version v2.4.1 of the code.
## PhaseTracer-1.0.2 [April 11, 2020]
 * Specify axis limits for phase_plotter
 * Update FS version used to 2.4.2, avoids compilation complaining about an unused Mathematica interafce problem on Mac OS with Mathamatica 12.
## PhaseTracer-1.0.3 [April 15, 2020]
 * Many thanks to Jingwei Lian for pointing out a bug in THDMIISNMSSMBCsimple.hpp. get_vector_debye_sq() returns two W boson masses, instead of one W boson mass and one Z boson mass.
## PhaseTracer-1.1.0 [January 4, 2021]
 * Update BSMPT used in PhaseTracer to version 2
 * Add requirement on cmake version, >=3.9, because of OpenMP support
 * Add a setting about BOOST in cmake to fix compiling problem with BOOST.1.72
 * Add functions, 'get_minima_at_t_low' and 'get_minima_at_t_high', to get minima at lowest and highest temperature. 
 * Add function, 'get_deepest_phase_at_T', to get the deepest phase at T.
 * Change 'V1' function in 'one_loop_potential' class to virtual function
 * Fix a bug that counter_term is not added to zero-temperature potential
## PhaseTracer-2.0 [November 5, 2021]
 * Added xi to OneLoopPotential class, as well as the relevant contributions to V1
 * Add high-temperature expansions into the one-loop potential class
 * Add On-shell like scheme example
 * Add covariant gauge example
 * Add h_bar_expansion method to obtain TC
 * Merge interface with TransitionSlover
 * Modify the one-loop potential class to interface with DRalgo
 * Add calculation of action, and relevant outputs
 * Add calculation of nucleation temperature
 * Add calculation of alpha, beta
 * Add calculation of GW spectrum, and relevant outputs
 * Add calculation of GW SNR, and relevant outputs
   

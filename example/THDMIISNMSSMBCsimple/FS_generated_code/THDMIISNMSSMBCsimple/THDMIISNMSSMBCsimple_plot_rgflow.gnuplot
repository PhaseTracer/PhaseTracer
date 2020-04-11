set title "THDMIISNMSSMBCsimple renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='THDMIISNMSSMBCsimple_rgflow.dat'

plot for [i=2:46+1] filename using 1:(column(i)) title columnhead(i)

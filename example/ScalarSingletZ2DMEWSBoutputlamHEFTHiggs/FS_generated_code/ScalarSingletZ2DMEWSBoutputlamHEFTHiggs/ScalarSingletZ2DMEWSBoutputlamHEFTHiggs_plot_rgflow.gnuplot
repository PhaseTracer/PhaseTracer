set title "ScalarSingletZ2DMEWSBoutputlamHEFTHiggs renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='ScalarSingletZ2DMEWSBoutputlamHEFTHiggs_rgflow.dat'

plot for [i=2:36+1] filename using 1:(column(i)) title columnhead(i)

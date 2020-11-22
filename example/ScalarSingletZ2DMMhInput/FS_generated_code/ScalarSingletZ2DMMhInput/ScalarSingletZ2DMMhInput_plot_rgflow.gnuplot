set title "ScalarSingletZ2DMMhInput renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='ScalarSingletZ2DMMhInput_rgflow.dat'

plot for [i=2:36+1] filename using 1:(column(i)) title columnhead(i)

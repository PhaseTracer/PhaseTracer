set title "ScalarSingletZ2DMMhInput particle spectrum"
set ylabel "mass / GeV"
unset key
unset bars

if (!exists("filename")) filename='ScalarSingletZ2DMMhInput_spectrum.dat'

plot filename using 1:2:(0.4):xtic(3) with xerrorbars pointtype 0 linecolor rgb "black"

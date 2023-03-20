set title "ScalarSingletZ2DMMhInputMsInput particle spectrum"
set ylabel "mass / GeV"
unset key
unset bars

if (!exists("filename")) filename='ScalarSingletZ2DMMhInputMsInput_spectrum.dat'

plot filename using 1:2:(0.4):xtic(3) with xerrorbars pointtype 0 linecolor rgb "black"

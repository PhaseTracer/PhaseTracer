set terminal pngcairo size 1024, 712 enhanced font 'Verdana,20'
set output png_name

set xlabel "x"
set ylabel  sprintf("V({/Symbol f}_f + x ({/Symbol f}_t - {/Symbol f}_f), T = %f)", T)

unset key
set autoscale xy

set xtics ('-0.25' -0.25, '0, {/Symbol f}_f' 0., '0.25' 0.25, '0.5' 0.5, '0.75' 0.75, '1, {/Symbol f}_t' 1., '1.25' 1.25)

plot dat_name u 1:2 w l

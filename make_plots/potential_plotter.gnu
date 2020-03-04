set terminal pngcairo size 1024, 712 enhanced font 'Verdana,20'
set output png_name

stats dat_name

set xlabel "{/Symbol f}_1"
set ylabel "{/Symbol f}_2"
set cblabel sprintf("V({/Symbol f}, T = %f)", T) 

unset key
set autoscale xy
set rmargin 10.

set xrange [STATS_min_x:STATS_max_x]
set yrange [STATS_min_y:STATS_max_y]

plot dat_name u 1:2:3 w image

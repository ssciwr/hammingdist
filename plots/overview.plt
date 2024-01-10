# to generate: gnuplot overview.plt
set terminal png enhanced crop
set output 'overview.png'
set title 'Hammingdist performance history (million genomes / second)'
set logscale y
set ylabel "million genomes / second"
set xrange [-0.3:4.3]
set yrange [0.005:200]
set bmargin 5
p 'overview.txt' u (1e3/$2):xtic(1) w lp t "" ps 2 pt 7 lc 1 lt 0

# to generate: gnuplot speed.plt
set terminal png enhanced crop
set output 'speed.png'
set title 'Hammingdist lower triangular output speed'
set xrange [0:330000]
set yrange [0:100]
set ylabel "Million distance matrix coefficients per second"
set xlabel "Sample count"
a = 0.5 * 1e-6
p \
'speed.txt' u 1:(a*($1**2)/$4) w lp t "from-fasta-to-lower-triangular GPU (A100)", \
'speed.txt' u 1:(a*($1**2)/$3) w lp t "from-fasta/dump-lower-triangular GPU (A100)", \
'speed.txt' u 1:(a*($1**2)/$2) w lp t "from-fasta/dump-lower-triangular CPU (52-core Xeon)"

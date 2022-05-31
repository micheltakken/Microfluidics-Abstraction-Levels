set datafile separator ','
set terminal pdf
set output "MixingProfiles.pdf"
set xlabel "Position across channel width in {/Symbol m}m"
set xrange [0:100]
set yrange [0:1]
set ylabel "Concentration of fluid A"
set style line 1 lc rgb '#0000C8' lt 1 lw 2 pt 4 ps 1.00 pi 100
set style line 2 lc rgb '#009600' lt 1 lw 2 pt 8 ps 1.00 pi 100
set style line 3 lc rgb '#C80000' lt 1 lw 2 pt 2 ps 1.00 pi 100
plot \
'1D/Fluid-Mixing/1D-Mixing.csv' using ($0*0.1):1 w lp ls 1 t "1D", \
'FVM/Fluid-Mixing/2D/Mixing2DFVM.csv' using ($0*0.1):4 w lp ls 2 t "2D FVM", \
'LBM/Fluid-Mixing/2D/Mixing2DLBM.csv'using ($0*0.1):1 w lp ls 3 t "2D LBM"

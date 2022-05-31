set datafile separator ','
set terminal pdf
set output "CarreauPlot.pdf"
set xlabel "Position across channel width in {/Symbol m}m"
set xrange [0:100]
set yrange [0:25]
set ylabel "Flow velocity in mm/s"
set style line 1 lc rgb '#000000' lt 1 lw 1 pt 4 ps 0.50 pi 80
set style line 2 lc rgb '#009600' lt 1 lw 1 pt 8 ps 0.70 pi 80
set style line 3 lc rgb '#0000C8' lt 1 lw 1 pt 2 ps 0.50 pi 100
set style line 4 lc rgb '#9400D3' lt 1 lw 1 pt 6 ps 0.50 pi 80
set style line 5 lc rgb '#C80000' lt 1 lw 1 pt 10 ps 0.70 pi 100
plot \
'1D/Non-Newtonian/1D-NonNewtonian.csv' using ($0*0.1):($1*1000) w lp ls 1 t "1D", \
'FVM/Non-Newtonian/2D/NonNewtonian-FVM2D.csv' using ($0*0.1):($1*1000) w lp ls 2 t "2D FVM", \
'FVM/Non-Newtonian/3D/NonNewtonian-FVM3D.csv' using ($0*0.1):($1*1000) w lp ls 3 t "3D FVM", \
'LBM/Non-Newtonian/2D/NonNewtonian-LBM2D.csv' using ($0*0.1):($2*1000) w lp ls 4 t "2D LBM", \
'LBM/Non-Newtonian/3D/NonNewtonian-LBM3D.csv' using ($0*0.1):($2*1000) w lp ls 5 t "3D LBM"


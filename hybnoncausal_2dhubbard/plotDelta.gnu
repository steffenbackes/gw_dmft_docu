#
set term post enhanced color "Helvetica" 22
set out "2dHubbard_Delta.ps"

xsize = 0.9
ysize = 0.9

set key samplen 2 at graph 0.99, graph 0.15
set xzeroaxis lt 1 lc 7 lw 1
set xlabel "{/Symbol t} [1/eV]"
set size xsize,2*ysize
set multiplot
set size xsize,ysize

set origin 0,ysize
set ylabel "{/Symbol D}({/Symbol t}) [eV]" offset 1
set label 1 "2D Hubbard model" at graph 0.12, graph 0.85
set label 2 "U=2, t=1, n=1" at graph 0.12, graph 0.75
set label 3 "{/Symbol D}({/Symbol t})" at graph 0.7, graph 0.9
p \
'delta_tauU2.dat'        u 1:2 w l lt 1 lw 3 lc 1 ti "non-causal" ,\
'delta_tauU2_causal.dat' u 1:2 w l lt 1 lw 3 lc 3 ti "causal" 

set origin 0,0
set label 2 "U=4, t=1, n=1" at graph 0.12, graph 0.75
p \
'delta_tau.dat'        u 1:2 w l lt 1 lw 3 lc 1 ti "non-causal" ,\
'delta_tau_causal.dat' u 1:2 w l lt 1 lw 3 lc 3 ti "causal" 

##
set size 0.5*xsize, 0.5*ysize
unset label 1
unset label 2
unset label 3
unset xlabel
unset ylabel
unset key
set xrange [39:40]
set xtics 0.5
set ytics 0.5

set origin 0.1*xsize, 1.15*ysize
p \
'delta_tauU2.dat'        u 1:2 w l lt 1 lw 3 lc 1 ti "non-causal" ,\
'delta_tauU2_causal.dat' u 1:2 w l lt 1 lw 3 lc 3 ti "causal" 

set origin 0.1*xsize, 0.15*ysize
p \
'delta_tau.dat'        u 1:2 w l lt 1 lw 3 lc 1 ti "non-causal" ,\
'delta_tau_causal.dat' u 1:2 w l lt 1 lw 3 lc 3 ti "causal" 

#
set term post enhanced color "Helvetica" 22
set out "Hybtest.ps"
set key samplen 2
set xlabel "{/Symbol w} [eV]"
set ylabel "{/Symbol D}({/Symbol w}) [eV]" offset 1
set xzeroaxis lt 1 lc 7

xsize = 0.9
ysize = 0.65

set size 1*xsize,3*ysize
set multiplot
set size xsize,ysize

set origin 0*xsize, 2*ysize
set label 1 "50% neg. weight!" at graph 0.05, graph 0.85
p \
'gloc.dat' u 1:(-$3/3.14) w l lt 1 lw 3 lc 1 ti 'A({/Symbol w})',\
'hybrid.dat' u 1:(-$3/3.14) w l lt 1 lw 3 lc 3 ti '{/Symbol D}({/Symbol w})'

set origin 0*xsize, 1*ysize
set key samplen 2 at graph 0.95, graph 0.7
set ytics 0.01
unset label 1
set xlabel "{/Symbol t} [1/eV]"
set ylabel "{/Symbol D}({/Symbol t}) [eV]" offset 1
p \
'hyb_tau.dat' u 1:2 w l lt 1 lw 3 lc 1 ti '{/Symbol D}({/Symbol t}) non-causal', \
'hyb_tau.dat' u 1:3 w l lt 1 lw 3 lc 3 ti '{/Symbol D}({/Symbol t}) causal'


set origin 0*xsize, 0*ysize
set ytics 0.5
set xrange [0:20]
set xlabel "i{/Symbol w}_n [eV]"
set ylabel "{/Symbol D}(i{/Symbol w}_n) [eV]" offset 1
p \
'hyb_mats.dat' u 1:3 w l lt 1 lw 3 lc 1 ti '{/Symbol D}(i{/Symbol w}_n) non-causal', \
'hyb_mats.dat' u 1:5 w l lt 1 lw 3 lc 3 ti '{/Symbol D}(i{/Symbol w}_n) causal'

#
unset key 
unset xlabel
unset ylabel
set size 0.5*xsize, 0.5*ysize

set origin 0.15*xsize, 1.2*ysize
set xrange [37:40]
set yrange [-0.05:0]
set xtics 1
set ytics 0.01
p \
'hyb_tau.dat' u 1:2 w l lt 1 lw 3 lc 1 ti '{/Symbol D}({/Symbol t}) non-causal', \
'hyb_tau.dat' u 1:3 w l lt 1 lw 3 lc 3 ti '{/Symbol D}({/Symbol t}) causal'


#
set term post enhanced color "Helvetica" 22

set out "Sigma_comp.ps"
set xlabel "i{/Symbol w}_n [eV]"
set ylabel "Im {/Symbol S}(i{/Symbol w}_n) [eV]" offset 1.5
set key samplen 2 at graph 0.95, graph 0.35

xsize = 0.75
ysize = 0.9
set size 2*xsize,1*ysize
set multiplot
set size xsize,ysize

set origin 0,0
set xrange [0:400]
set label 2 "SrVO_3" at graph 0.7, graph 0.9
p \
'sigma_GW.dat'  u 1:3 w l lt 1 lw 4 lc 1 ti "GW_{loc}" , \
'sigma_loc.dat' u 1:3 w l lt 1 lw 4 lc 2 ti "GW (cRPA screened)" , \
'sigma_imp.dat' u 1:3 w l lt 1 lw 4 lc 3 ti "GW (imp. screened)"

set origin xsize,0
set key samplen 2 at graph 0.75, graph 0.25
set xrange[0:2]
p \
'sigma_GW.dat'  u 1:3 w l lt 1 lw 4 lc 1 ti "GW_{loc}" , \
'sigma_loc.dat' u 1:3 w l lt 1 lw 4 lc 2 ti "GW (cRPA screened)" , \
'sigma_imp.dat' u 1:3 w l lt 1 lw 4 lc 3 ti "GW (imp. screened)"

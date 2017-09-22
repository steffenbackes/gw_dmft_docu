#
set term post enhanced color "Helvetica" 22

set out "W_loc_comp.ps"
set xrange [0:14]
set xlabel "i{/Symbol n}_n [eV]"
set ylabel "W(i{/Symbol n}_n) [eV]" offset 1.5
set key samplen 2 at graph 0.95, graph 0.35

xsize = 0.75
ysize = 0.9
set size 2*xsize,1*ysize
set multiplot
set size xsize,ysize

set origin 0,0
set label 1 "W_{abab} = U_0^{scr}" at graph 0.1, graph 0.8
set label 2 "SrVO_3" at graph 0.1, graph 0.9
p \
'Jmatrix_bare_crpa.dat' u 1:2 w l lt 1 lw 4 lc 1 ti "cRPA (bare)" ,\
'Jmatrix_scr_GW.dat'  u 1:2 w l lt 1 lw 4 lc 2 ti "W_{loc}" , \
'Jmatrix_bare_GW.dat' u 1:2 w l lt 1 lw 4 lc 3 ti "W_{loc} unscreened" ,\
'Jmatrix_scr_loc.dat' u 1:2 w l lt 1 lw 4 lc 4 ti "cRPA screened" , \
'Jmatrix_scr_imp.dat' u 1:2 w l lt 1 lw 4 lc 5 ti "U_{imp} screened" 

set origin xsize,0
set label 1 "W_{abba} = J^{scr}" at graph 0.65, graph 0.65
p \
'Jmatrix_bare_crpa.dat' u 1:3 w l lt 1 lw 4 lc 1 ti "cRPA (bare)" ,\
'Jmatrix_scr_GW.dat'  u 1:3 w l lt 1 lw 4 lc 2 ti "W_{loc}" , \
'Jmatrix_scr_loc.dat' u 1:3 w l lt 1 lw 4 lc 3 ti "cRPA screened" , \
'Jmatrix_scr_imp.dat' u 1:3 w l lt 1 lw 4 lc 4 ti "U_{imp} screened" 

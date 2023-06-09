#
set term post enhanced color "Helvetica" 22
set out "2dHubbard_Hyb.ps"

xsize = 0.9
ysize = 0.9

set key samplen 2
set xzeroaxis lt 1 lc 7 lw 1
set xlabel "{/Symbol w} [eV]"
set size 2*xsize,1*ysize
set multiplot
set size xsize,ysize

set origin 0,0*ysize
set ylabel "-Im G_{loc}({/Symbol w}) [eV]" offset 1
set label 1 "2D Hubbard model" at graph 0.02, graph 0.95
set label 2 "U=4, t=1, n=1" at graph 0.02, graph 0.85
set label 3 "A({/Symbol w})" at graph 0.7, graph 0.55
set yrange [0:0.25]
set mytics 2
p \
'G_loc.dat' u 1:(-$3/3.141) w l lt 1 lw 3 lc 1 ti "{/Symbol S}^{(2)}+DMFT" ,\
'G_locdmft.dat'  u 1:(-$3/3.141) w l lt 1 lw 3 lc 3 ti "DMFT" 


set origin xsize,0
set ylabel "Im. Hybridization {/Symbol D}({/Symbol w}) [eV]" offset 1
set label 1 "2D Hubbard model" at graph 0.02, graph 0.95
set label 2 "U=4, t=1, n=1" at graph 0.02, graph 0.85
set label 3 "{/Symbol D}({/Symbol w})" at graph 0.7, graph 0.55
set yrange [-0.05:0.65]
set mytics 2
p \
'Hyb_fullk.dat' u 1:(-$3/3.141) w l lt 1 lw 3 lc 1 ti "{/Symbol S}^{(2)}+DMFT" ,\
'Hyb_dmft.dat'  u 1:(-$3/3.141) w l lt 1 lw 3 lc 3 ti "DMFT" 

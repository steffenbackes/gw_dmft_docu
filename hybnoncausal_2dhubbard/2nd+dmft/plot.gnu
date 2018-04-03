#
set term post enhanced color "Helvetica" 22
set out "dmft+2nd.ps"

xsize = 0.9
ysize = 0.9

set key samplen 2 at graph 1, graph 0.4
set xzeroaxis lt 1 lc 7 lw 1
set xlabel "i{/Symbol w}_n [eV]"
set xrange [0:20]

set size xsize,2*ysize
set multiplot
set size xsize,ysize

set origin 0,ysize
set ylabel "Re {/Symbol S}_{loc}(i{/Symbol w}_n) [eV]" offset 1
set label 1 "2D Hubbard model" at graph 0.12, graph 0.95
set label 2 "U=4, t=1, n=0.8" at graph 0.12, graph 0.85
p \
'diagMCreal.dat'     u 1:2              w p pt 7 ps 1.5 lc 7 ti "diag.MC" ,\
'sigma_2nd.dat'      u 1:($2-4*0.4)     w l lt 1 lw 3 lc 1 ti "{/Symbol S}^{(2)}", \
'sigma_2ndLOCAL.dat' u 1:(($2-1.605)*2) w l lt 1 lw 3 lc 2 ti "loc {/Symbol S}^{(2)}", \
'sigma_dmft.dat'     u 1:($2-4*0.4)     w l lt 1 lw 3 lc 3 ti "DMFT", \
'sigma_2nd+dmft.dat' u 1:($2-4*0.4)     w l lt 1 lw 3 lc 4 ti "{/Symbol S}^{(2)}(0)+DMFT"

set origin 0,0
set ylabel "Im {/Symbol S}_{loc}(i{/Symbol w}_n) [eV]" offset 1
set label 1 "2D Hubbard model" at graph 0.12, graph 0.95
set label 2 "U=4, t=1, n=0.8" at graph 0.12, graph 0.85
p \
'diagMCimag.dat'     u 1:2      w p pt 7 ps 1.5 lc 7 ti "diag.MC" ,\
'sigma_2nd.dat'      u 1:3      w l lt 1 lw 3 lc 1 ti "{/Symbol S}^{(2)}", \
'sigma_2ndLOCAL.dat' u 1:($3*2) w l lt 1 lw 3 lc 2 ti "loc. {/Symbol S}^{(2)}", \
'sigma_dmft.dat'     u 1:3      w l lt 1 lw 3 lc 3 ti "DMFT", \
'sigma_2nd+dmft.dat' u 1:3      w l lt 1 lw 3 lc 4 ti "{/Symbol S}^{(2)}(0)+DMFT"

#
set term post enhanced color "Helvetica" 22
set size 0.85

set key samplen 2

set xlabel "i{/Symbol w} [eV]"
set label 1 "U=2" at graph 0.05, graph 0.95
set out "1dhub_simp.ps"

set ylabel "Im {/Symbol S}(i{/Symbol w}) [eV]"
set xrange [0:10]
p \
'sigk0/Swl.dat'    u 1:3 w l lt 1 lw 3 lc 7 ti "nonint", \
'gen_sigk/Swl.dat' u 1:3 w l lt 1 lw 3 lc 1 ti "Fock general", \
'dys_sigk/Swl.dat' u 1:3 w l lt 1 lw 3 lc 3 ti "Fock dyson"


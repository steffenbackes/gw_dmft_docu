#
set term post enhanced color "Helvetica" 22
set size 0.85

set key samplen 2

set xlabel "{/Symbol w} [eV]"

set out "1dhub_bath_sigk.ps"

set ylabel "Gbath({/Symbol w})"
set xrange [-6:6]
p \
'bathhybK.dat'         u 1:6 w l lt 1 lw 3 lc 7 ti "nonint", \
'bathhybKmultiorb.dat' u 1:2 w l lt 1 lw 3 lc 1 ti "Fock general", \
''                     u 1:4 w l lt 1 lw 3 lc 3 ti "Fock dyson"

set out "1dhub_hyb_sigk.ps"

set ylabel "{/Symbol D}({/Symbol w})"
set xrange [-6:6]
p \
'bathhybK.dat'         u 1:7 w l lt 1 lw 3 lc 7 ti "nonint", \
'bathhybKmultiorb.dat' u 1:3 w l lt 1 lw 3 lc 1 ti "Fock general", \
''                     u 1:5 w l lt 1 lw 3 lc 3 ti "Fock dyson"

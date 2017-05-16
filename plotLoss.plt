set term eps
set output "constComps.eps"


plot "Q_cGTY.dat" w l lw 2 lc rgb "blue" smooth csplines, \
 "Q_cT.dat" w l lw 2 lc rgb "red" smooth csplines, \
 "Q_cN.dat" w l lw 2 lc rgb "purple" smooth csplines, \
 "Exp.dat" u 1:(1e3 * $2) w l lw 2 lc rgb "black", \
 "Bulk.dat" u 1:(0.1 * $2) w l lt 10 lw 2 lc rgb "black"

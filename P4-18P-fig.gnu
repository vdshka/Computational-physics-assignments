
# Figure 1- Isotherm which contains only stable points (V,P) (Condition P'(V)<0)

set term png
set output "P4-18P-fig1.png"
set key top right
set xlabel "v(reduced units)"
set ylabel "P,dP(reduced units)"
plot "P4-18P-res.dat" index 0 u 1:($2<0?$3:1/0) w l title "P(v)"


# Figure 2- Curve at T=0.92 with V€[1/3+0.1,2]

set term png
set output "P4-18P-fig2.png"
set key top right
set xlabel "v(reduced units)"
set ylabel "P(reduced units)"
plot "P4-18P-res.dat" index 1 u 1:2 w l lc rgb "#CF3476" title "P(v)",\
0 notitle lc black
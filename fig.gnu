

# Figure 1- Sum depending on N and asymptotic behaviour


set term png
set output "P1-18P-fig1.png"
set key bottom right
set xlabel "N"
set ylabel "SUM(8,N)"
set logscale y
plot "P1-18P-res1.dat" u 1:2 pt 7 title "SUM(8,N) (N)",\
"P1-18P-res1.dat" u 1:2 w l title "asymptotic behaviour"



# Figure 2- Sum divided by the asymptotic behaviour depending on N 


set term png
set output "P1-18P-fig2.png"
set key top right
set xlabel "N"
set ylabel "SUM(8,N)/SUM asympt"
plot "P1-18P-res2.dat" u 1:2 pt 7 title "SUM(8,N)/SUM asympt (N)"

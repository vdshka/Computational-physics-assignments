
# Figure 1- Error depending on the length of the intervals. Comparison with the expected behaviour.

set term png
set output "P3-18P-fig1.png"
set key top left
set logscale y
set logscale x
set format y "10^{%+L}"
set format x "10^{%+L}"
set xlabel "h (10^6 km)"
set ylabel "error(10^1^2 km^2)"
plot "P3-18P-res1.dat" u 1:4 pt 7 title " Error trapezium",\
"P3-18P-res1.dat" u 1:5 pt 7 title " Error Simpson",\
"P3-18P-res1.dat" u 1:6 w l title "Expected trapezium",\
"P3-18P-res1.dat" u 1:7 w l title "Expected Simpson"



# Figure 2- Convergence compared with the expected behaviour.

set term png
set output "P3-18P-fig2.png"
set key top left
set logscale y
set logscale x
set format y "10^{%+L}"
set format x "10^{%+L}"
set xlabel "h (10^6 km)"
set ylabel "error(10^1^2 km^2)"
plot "P3-18P-res2.dat" u 1:4 pt 7 title " Error trapezium",\
"P3-18P-res2.dat" u 1:5 pt 7 title " Error Simpson",\
"P3-18P-res2.dat" u 1:6 w l title "Expected trapezium",\
"P3-18P-res2.dat" u 1:7 w l title "Expected Simpson"
#----------------------------------------------------------------------------------

# Figure 1- Positions of the first and fifth piston depending on time

#----------------------------------------------------------------------------------
set term png
set output "P2-18P-fig1.png"
set key top left
set xlabel "time(s)"
set ylabel "piston position(cm)"
plot "P2-18P-res1.dat" u 1:2 pt 7 title " P1",\
"P2-18P-res1.dat" u 1:6 pt 7 title " P5"

#-----------------------------------------------------------------------------------------------

# Figure 2- Positions of the third and fifth piston depending on the position of the second one

#-----------------------------------------------------------------------------------------------
set term png
set output "P2-18P-fig2.png"
set key top left
set xlabel "Position second piston(cm)"
set ylabel "Position piston(cm)"
plot "P2-18P-res1.dat" u 3:4 pt 7 title " P3",\
"P2-18P-res1.dat" u 3:6 pt 7 title " P5"

#-------------------------------------------------------------------------------------------------------

# Figure 3- Position of the fourth piston computed directly and with interpolation (order zero and linear) 

#-------------------------------------------------------------------------------------------------------
set term png
set output "P2-18P-fig3.png"
set key top right
set xlabel "t (s)"
set ylabel "Fourth piston position(cm)"
plot "P2-18P-res2.dat" u 1:2 w l title "int zero",\
"P2-18P-res2.dat" u 1:3 pt 3 ps 0 title "int lin",\
"P2-18P-res1.dat" u 1:5 pt 7 ps 1 title "X4 (t)"
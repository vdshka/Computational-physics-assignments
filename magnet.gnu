#----------------------------------------------------------------------------------

# Figura 1- magnetització en funció del temps temperatura= 2

#----------------------------------------------------------------------------------
set term png
set output "fig1.png"
set key bottom right box 
set title "Temperatura=2"
set xlabel "Magnetització"
set ylabel "Temps"
set grid
plot "magnettemp1.dat" u 1:2 w l title "Magnetització(temps)"

#----------------------------------------------------------------------------------

# Figura 2- magnetització en funció del temps temperatura= 3.5

#----------------------------------------------------------------------------------
set term png
set output "fig2.png"
set key bottom right box 
set title "Temperatura=3.5"
set xlabel "Magnetització"
set ylabel "Temps"
set grid
plot "magnettemp2.dat" u 1:2 w l title "Magnetització(temps)"

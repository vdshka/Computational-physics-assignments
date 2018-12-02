#----------------------------------------------------------------------------------

# Figura 1- Grafica 2D amb fogons

#----------------------------------------------------------------------------------

set xlabel "x (cm)"
set ylabel "y (cm)"
set palette rgbformulae 30,31,32
set title "Mapa de temperatures amb fogons"
set term png
set view map
set size ratio -1
set output "P9-18P-fig4.png"

splot "Amb_fogons.dat" u 1:2:3 w pm3d t""


#----------------------------------------------------------------------------------

# Figura 2- Grafica 2D sense fogons

#----------------------------------------------------------------------------------

set xlabel "x (cm)"
set ylabel "y (cm)"
set palette rgbformulae 30,31,32
set title "Mapa de temperatures sense fogons"
set term png
set size ratio -1
set view map
set output "P9-18P-fig5.png"

splot "Sense_fogons.dat" u 1:2:3 w pm3d t""

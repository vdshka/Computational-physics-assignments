#----------------------------------------------------------------------------------

# Figura 1- Comparacio dels dos metodes per a Tinterior = 10ºC

#----------------------------------------------------------------------------------

set term png
set output "P9-18P-fig1.png"
set key bottom right box 
set title "Tinterior = 10 ºC"
set xlabel "Nº iteracions"
set ylabel "Temperatura (ºC)"
set grid
plot "P9-18P.dat" index 0 u 1:2 w l title "Mètode Gauss-Seidel",\
"P9-18P.dat" index 1 u 1:2 w l title "Mètode de Sobrerelaxació"

#----------------------------------------------------------------------------------

# Figura 2- Comparacio dels dos metodes per a Tinterior = 120ºC

#----------------------------------------------------------------------------------

set output "P9-18P-fig2.png"
set key top right box 
set title "Tinterior = 120 ºC"
set xlabel "Nº iteracions"
set ylabel "Temperatura (ºC)"
set grid
plot "P9-18P.dat" index 2 u 1:2 w l title "Mètode Gauss-Seidel",\
"P9-18P.dat" index 3 u 1:2 w l title "Mètode de Sobrerelaxació"

#----------------------------------------------------------------------------------

# Figura 3- Comparacio dels dos metodes per a Tinterior = 1040ºC

#----------------------------------------------------------------------------------

set output "P9-18P-fig3.png"
set key top right box 
set title "Tinterior = 1040 ºC"
set xlabel "Nº iteracions"
set ylabel "Temperatura (ºC)"
set grid
plot "P9-18P.dat" index 4 u 1:2 w l title "Mètode Gauss-Seidel",\
"P9-18P.dat" index 5 u 1:2 w l title "Mètode de Sobrerelaxació"


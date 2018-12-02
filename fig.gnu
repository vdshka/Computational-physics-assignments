#----------------------------------------------------------------------------------

# Figura 1- Solucions en interval x € [-L/2:L/4]

#----------------------------------------------------------------------------------
set term png
set output "P8-18P-fig1.png"
set key top left box 
set title "500 passos "
set xlabel "x"
set xrange [-3.5:1.75]
set ylabel "phi(x,E)"
set grid
plot "P8-18P-res1.dat" index 0 u 1:2 w p pt 2 title "E1",\
"P8-18P-res1.dat" index 1 u 1:2 w p pt 2 title "E2",\
"P8-18P-res1.dat" index 2 u 1:2 w p pt 2 title "E3",\
"P8-18P-res1.dat" index 3 u 1:2 w p pt 2 title "E4"
unset xrange

#----------------------------------------------------------------------------------

# Figura 2- Convergència. Energia a cada iteració pels 3 autovalors.

#----------------------------------------------------------------------------------
set term png
set output "P8-18P-fig2.png"
set key top left box 
set title "500 passos "
set xlabel "x"
set ylabel "phi(x,E)"
set grid
plot "Fig.dat" index 0 u 1:2 w p pt 2 title "E1",\
"Fig.dat" index 1 u 1:2 w p pt 2 title "E2",\
"Fig.dat" index 2 u 1:2 w p pt 2 title "E3"

#----------------------------------------------------------------------------------

# Figura 3- Autovectors normalitzats.

#----------------------------------------------------------------------------------
set term png
set output "P8-18P-fig3.png"
set key top right box 
set title "500 passos "
set xlabel "Energia"
set ylabel "phi(x=L/2)"
set grid
plot "P8-18P-res3.dat" index 0 u 1:2 w l title "Autovalors",\
"P8-18P-res3.dat" index 0 u 1:2 w p notitle

#----------------------------------------------------------------------------------

# Figura 4- Estat fonamental normalitzat per a diferents valors de beta.

#----------------------------------------------------------------------------------
set term png
set output "Fig4.png"
set key top left box 
set xlabel "x"
set ylabel "phi(x,E)"
set grid
plot "Fig.dat" index 0 u 1:2 w l title "Beta=0",\
"Fig.dat" index 3 u 1:2 w l lc rgb "orange"  title "Beta=1",\
"Fig.dat" index 6 u 1:2 w l title "Beta=5"

#----------------------------------------------------------------------------------

# Figura 5- Probabilitat

#----------------------------------------------------------------------------------
set term png
set output "Fig5.png"
set key top right box 
set xlabel "x"
set ylabel "phi(x,E)"
set xrange [-1.5:2]
set grid
plot "Probabilitat.dat" index 0 u 1:2 w l title "Beta=0,E=-79.09",\
"Probabilitat.dat" index 1 u 1:2 w l title "Beta=0,E=-76.36",\
"Probabilitat.dat" index 2 u 1:2 w l title "Beta=0,E=-71.81",\
"Probabilitat.dat" index 3 u 1:2 w l title "Beta=1,E=-78.00",\
"Probabilitat.dat" index 4 u 1:2 w l title "Beta=1,E=-73.73",\
"Probabilitat.dat" index 5 u 1:2 w l title "Beta=1,E=68.60",\
"Probabilitat.dat" index 6 u 1:2 w l title "Beta=5,E=-75.64",\
"Probabilitat.dat" index 7 u 1:2 w l title "Beta=5,E=-75.64",\
"Probabilitat.dat" index 8 u 1:2 w l title "Beta=5,E=-66.90"

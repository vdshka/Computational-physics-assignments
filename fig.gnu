#----------------------------------------------------------------------------------

# Figura 1- perfil estacionari

#----------------------------------------------------------------------------------
set term png
set output "P10-18P-fig1.png"
set key bottom right box 
set title "Situació inicial"
set xlabel "x (cm)"
set ylabel "Temperatura (ºC)"
set grid
plot "P10-18P.dat" index 0 u 1:2 w l title "Gauss-Seidel",\
"P10-18P.dat" index 1 u 1:2 w l title "Jacobi",\
"P10-18P.dat" index 2 u 1:2 w l title "Sobrerelaxació"

#-----------------------------------------------------------------------------------

# Figura 2- Evolució temporal temperatures xp= 6,10,28

#-----------------------------------------------------------------------------------
set term png
set output "P10-18P-fig2.png"
set key top right box  
set title "Evolució T(xp,t)(ºC)"
set xlabel "t (s)"
set ylabel "Temperatura (ºC)"
set yrange [15:55]
set grid
plot "P10-18P.dat" index 3 u 1:2 w l title "xp =  6  cm",\
"P10-18P.dat" index 3 u 1:3 w l title "xp = 10 cm",\
"P10-18P.dat" index 3 u 1:4 w l title "xp = 16 cm",\
"P10-18P.dat" index 3 u 1:5 w l title "xp = 28 cm"

unset yrange
#----------------------------------------------------------------------------------- 

# Figura 3- Evolucio temporal de la temperatura mitjana

#-----------------------------------------------------------------------------------
set term png
set output "P10-18P-fig3.png"
set key top right box 
set title "Evolució Tmitjana(t)(ºC)"
set xlabel "t (s)"
set ylabel "Tmitjana (ºC)"
set grid
plot "P10-18P.dat" index 4 u 1:2 w l title "k =   3.0",\
"P10-18P.dat" index 5 u 1:2 w l title "k =   7.4",\
"P10-18P.dat" index 6 u 1:2 w l title "k = 10.0"

#-----------------------------------------------------------------------------------

# Figura 4- Evolucio temporal perfil de temperatures per a k=10
# per veure que la evolució té bona pinta ^^

#-----------------------------------------------------------------------------------
set term png
set output "P10-18P-fig4.png"
#datafile ="P10-18P.dat"
set title "Evolució temporal perfil de temperatures, k = 10"
set xlabel "x(cm)"
set ylabel "Temperatura (ºC)"

plot for [idx=7:125] "P10-18P.dat" i idx u 1:2 w l notitle

#------------------------------------------------------------------------------------

# Animació evolució temporal perfil de temperatures per a k=10

#------------------------------------------------------------------------------------

set term gif animate delay 20 
set output "anima.gif"
set title "Evolució temporal perfil de temperatures, k = 10"
set key top left box
set xrange[0:30]
set yrange[10:50]
set xlabel "x (cm)"
set ylabel "Temperatura (ºC)"

do for[i=7:125]{
plot "P10-18P.dat" index 7 u 1:2 w l title "Perfil inicial   ",\
"P10-18P.dat" index i u 1:2 w l title 'Temps = '.columnhead(2)." s"
}


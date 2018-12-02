#----------------------------------------------------------------------------------

# Figura 1- espins

#----------------------------------------------------------------------------------
set term png size 500,500
set output n.".png"
set size square
unset key
set xrange[0:51]
set yrange[0:51]
set xlabel "coordenada i"
set ylabel "coordenada j"
plot n.".dat" u 1:2 w p pt 5 ps 0.8
set output
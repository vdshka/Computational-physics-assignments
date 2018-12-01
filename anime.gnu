set term gif size 1000,600 animate  delay 20 loop 0 
set output "P7-ani-b.gif"
datafile ="anime.dat"


do for[i=1:1000:5]{

set multiplot

set size 0.5,0.8
set origin 0.0,0.0
set title "As a function of time"
set size 0.5,0.8
set xrange[0:50]
set yrange[-2*pi:2*pi]
set xlabel "t (s)"
set ylabel "theta, thetap"
plot datafile every ::1::i with line linewidth 4 t"theta" ,datafile every ::1::i u 1:4 with line linewidth 4 t"thetap"

set origin 0.5,0
set size 0.5,0.8
set title "Phase diagram"
set yrange[-2*pi:2*pi]
set xrange[-pi:pi]
set xlabel "theta"
set ylabel "thetap"
if (i>200) { plot datafile every::i::i u 2:3 t"" ps 3,datafile every::i-100::i u 4:5 w l t"" } else {
plot datafile every::i::i u 2:3 t"" ps 3}
unset multiplot

}

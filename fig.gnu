

# Figure 1- Small-angle motion.

set term png
set output "P7-18P-b-fig1.png"
set key bottom box outside horizontal
set title "Small-angle - angle vs t"
set xlabel "Time (s)"
set ylabel "Angle (rad)"
set grid
plot "P7-18P-res.dat" index 0 u 1:2 w l lc rgb "orange" title "Euler raw",\
"P7-18P-res.dat" index 0 u 1:4 w l title "Predictor/corrector",\
"P7-18P-res.dat" index 1 u 1:2 w l lc rgb "#CC0605" title "Approximation"


# Figure 2- Large-angle motion. Pendulum dynamics.

set term png
set output "P7-18P-b-fig2.png"
set key bottom left box
set title "Large-angle - angle vs t"
set xlabel "Time (s)"
set ylabel "Angle (rad)"
set grid
plot "P7-18P-res.dat" index 2 u 1:2 w l lc rgb "orange" title "Euler raw",\
"P7-18P-res.dat" index 2 u 1:4 w l title "Predictor/corrector"


# Figure 3- Large-angle motion. Trajectories in the phase space.

set term png
set output "P7-18P-b-fig3.png"
set key top left box
set title "Large-angle - trajectories"
set xlabel "Angle (rad)"
set ylabel "Velocity (rad/s)"
set grid
plot "P7-18P-res.dat" index 2 u 2:3 w l lc rgb "orange" title "Euler raw",\
"P7-18P-res.dat" index 2 u 4:5 w l title "Predictor/corrector"


# Figure 4- Evolution of the kinetic and total energy. 


set term png
set output "P7-18P-b-fig4.png"
set key top left box 
set title "KINETIC ENERGY (K) AND TOTAL ENERGY (E)"
set xlabel "Time (s)"
set ylabel "Energy (J)"
set grid
plot "P7-18P-res.dat" index 3 u 1:4 w l lc rgb "orange" title "E(t) euler raw ",\
"P7-18P-res.dat" index 3 u 1:7 w l title "E(t) predictor/corrector ",\
"P7-18P-res.dat" index 3 u 1:2 w l lc rgb "#EA899A" title "K(t) euler raw ",\
"P7-18P-res.dat" index 3 u 1:5 w l lc rgb "#8673A1" title "K(t) predictor/corrector "


# Figure 5- Transition. Trajectories in the phase space.


set term png
set output "P7-18P-b-fig5.png"
set key bottom right box
set title "TRANSITION - trajectories"
set xlabel "Angle (rad)"
set ylabel "Velocity (rad/s)"
set grid
plot "P7-18P-res.dat" index 4 u 2:3 w l lc rgb "#8673A1" title "alpha = 2sqrt(g/l)+0.1",\
"P7-18P-res.dat" index 4 u 4:5 w l title "alpha = 2sqrt(g/l)-0.1"


# Figure 6- Convergence of the method. Total energy as a function of time.


set term png
set output "P7-18P-b-fig6.png"
set key top left box
set title "CONVERGENCE - E(t)"
set xlabel "Time(t)"
set ylabel "Total energy (J)"
set grid
plot "P7-18P-res.dat" index 5 u 1:2 w l lc rgb "#49678D" title "200 passos de temps",\
"P7-18P-res.dat" index 6 u 1:2 w l lc rgb "red" title "600 passos de temps",\
"P7-18P-res.dat" index 7 u 1:2 w l lc rgb "green" title "4000 passos de temps",\
"P7-18P-res.dat" index 8 u 1:2 w l lc rgb "black" title "50000 passos de temps"
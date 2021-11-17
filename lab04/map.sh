set term png

set output "map_1.0.png"
set xlabel "x"
set ylabel "y"
set title "Mapa V(x,y) Ω=1.0"
set pm3d map
splot [0:15][0:10] "map_1.0.dat" i 0 u 1:2:3

set output "map_0.6.png"
set xlabel "x"
set ylabel "y"
set title "Mapa V(x,y) Ω=0.6"
set pm3d map
splot [0:15][0:10] "map_0.6.dat" i 0 u 1:2:3


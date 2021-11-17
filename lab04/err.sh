set terminal png large size 960,540
set output "error_1.0.png"
set xlabel "x"
set ylabel "y"
set size ratio -1
set title "Blad rozwiazania dla relaksacji globalnej Ω=1.0"
set pm3d map
set palette rgbformulae 33,13,10
splot [0:15][0:10] "error_1.0.dat" i 0 u 1:2:3


set output "error_0.6.png"
set xlabel "x"
set ylabel "y"
set size ratio -1
set title "Blad rozwiazania dla relaksacji globalnej Ω=0.6"
set pm3d map
set palette rgbformulae 33,13,10
splot [0:15][0:10] "error_0.6.dat" i 0 u 1:2:3


#!/usr/bin/gnuplot
set term png

set output "map_n_50.png"
set xlabel "x"
set ylabel "y"
set title "nx = ny = 50"
set pm3d map
set palette defined (-10 "blue", 0 "green", 10 "yellow")
set size ratio -1

splot [0:5][0:5] "map_n_50.txt" i 0 u 2:1:3

reset

set output "map_n_100.png"
set xlabel "x"
set ylabel "y"
set title "nx = ny = 100"
set pm3d map
set palette defined (-10 "blue", 0 "green", 10 "yellow")
set size ratio -1

splot [0:10][0:10] "map_n_100.txt" i 0 u 2:1:3

reset

set output "map_n_200.png"
set xlabel "x"
set ylabel "y"
set title "nx = ny = 200"
set pm3d map
set palette defined (-10 "blue", 0 "green", 10 "yellow")
set size ratio -1

splot [0:20][0:20] "map_n_200.txt" i 0 u 2:1:3

reset

set output "map_eps_1.png"
set xlabel "x"
set ylabel "y"
set title "epsilon_1 = epsilon_2 = 1"
set pm3d map
set palette defined (-10 "blue", 0 "green", 10 "yellow")
set size ratio -1

splot [0:10][0:10][-0.8:0.8] "map_eps_1.txt" i 0 u 2:1:3


reset

set output "map_eps_1_2.png"
set xlabel "x"
set ylabel "y"
set title "epsilon_1 = 1, epsilon_2 = 2"
set pm3d map
set palette defined (-10 "blue", 0 "green", 10 "yellow")
set size ratio -1
set cbrange [-0.8:0.8]

splot [0:10][0:10][-0.8:0.8] "map_eps_1_2.txt" i 0 u 2:1:3


reset

set output "map_eps_1_10.png"
set xlabel "x"
set ylabel "y"
set title "epsilon_1 = 1, epsilon_2 = 10"
set pm3d map
set palette defined (-10 "blue", 0 "green", 10 "yellow")
set size ratio -1
set cbrange [-0.8:0.8]

splot [0:10][0:10][-0.8:0.8] "map_eps_1_10.txt" i 0 u 2:1:3
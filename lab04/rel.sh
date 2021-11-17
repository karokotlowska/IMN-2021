set term png

set output "global.png"
set logscale x
set xlabel "it"
set ylabel "S(it)"
set grid
set title "Całka S(it) dla relaksacji globalnej"
plot "global_1.0.dat" u 1:2 w l t "Ω = 1.0" , "global_0.6.dat" u 1:2 w l lw 2 t "Ω = 0.6"

set output "local.png"
set logscale x
set logscale x
set xlabel "it"
set ylabel "S(it)"
set title "Całka S(it) dla relaksacji lokalnej"
plot "local_1.0.dat" u 1:2 w l t "Ω = 1.0" , "local_1.4.dat" u 1:2 w l lw 2 t "Ω = 1.4", "local_1.8.dat" u 1:2 w l lw 2 t "Ω = 1.8", "local_1.9.dat" u 1:2 w l lw 2 t "Ω = 1.9"

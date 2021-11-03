set term png 
set border linewidth 2
set xl "u(t)"
set yl "u(t),z(t)"
set linetype 1 linewidth 3
set linetype 2 linewidth 3

set out "picard.png"
set title "Metoda Picarda"
p "pickard.txt" u 1:2 w l t "u(t)", \
  "pickard.txt" u 1:3 w l t "z(t)=N-u(t)"

set out "newton.png"
set title "Iteracja Newtona"
p "newton.txt" u 1:2 w l t "u(t)", \
  "newton.txt" u 1:3 w l t "z(t)=N-u(t)"

set out "RK2.png"
set title "Niejawna RK2"
p "RK2.txt" u 1:2 w l t "u(t)", \
  "RK2.txt" u 1:3 w l t "z(t)=N-u(t)"

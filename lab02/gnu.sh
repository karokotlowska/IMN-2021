set term png 

set xl "t"
set yl "u(t), v(t)"

set out "metoda_picarda.png"
set title "Metoda Picarda"
p "pickard.txt" u 1:2 w l t "u(t)", \
  "pickard.txt" u 1:3 w l t "v(t)"

set out "iteracja_newtona.png"
set title "Iteracja Newtona"
p "newton.txt" u 1:2 w l t "u(t)", \
  "newton.txt" u 1:3 w l t "v(t)"

set out "rk_2.png"
set title "Niejawna RK2"
p "RK2.txt" u 1:2 w l t "u(t)", \
  "RK2.txt" u 1:3 w l t "v(t)"

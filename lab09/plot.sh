set term png size 500,500 


set xl "x"
set yl "y"

set view map



set title "T(it=100)"
set out "T_it=100.png"
splot 'T.dat' i 0 u 2:3:4 w pm3d notitle

set title "T(it=200)"
set out "T_it=200.png"
splot 'T.dat' i 1 u 2:3:4 w pm3d notitle

set title "T(it=500"
set out "T_it=500.png)"
splot 'T.dat' i 2 u 2:3:4 w pm3d notitle

set title "T(it=1000)"
set out "T_it=1000.png"
splot 'T.dat' i 3 u 2:3:4 w pm3d notitle

set title "T(it=2000)"
set out "T_it=2000.png" 
splot 'T.dat' i 4 u 2:3:4 w pm3d notitle




#set title "∇^2T(x,y) it=100"
set title "grad^2T(it=100)"
set out "pow_nabla_T_it=100.png"
splot 'pow_nabla_T.dat' i 0 u 2:3:4 w pm3d notitle

#set title "∇^2T(x,y) it=200"
set title "grad^2T(it=200)"
set out "pow_nabla_T_it=200.png"
splot 'pow_nabla_T.dat' i 1 u 2:3:4 w pm3d notitle

#set title "∇^2T(x,y) it=500"
set title "grad^2T(it=500)"
set out "pow_nabla_T_it=500.png"
splot 'pow_nabla_T.dat' i 2 u 2:3:4 w pm3d notitle

#set title "∇^2T(x,y) it=1000"
set title "grad^2T(it=1000)"
set out "pow_nabla_T_it=1000.png"
splot 'pow_nabla_T.dat' i 3 u 2:3:4 w pm3d notitle

#set title "∇^2T(x,y) it=2000"
set title "grad^2T(it=2000)"
set out "pow_nabla_T_it=2000.png" 
splot 'pow_nabla_T.dat' i 4 u 2:3:4 w pm3d notitle

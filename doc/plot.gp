set term x11
set macro
u = 'u 2:($3/75) w l lw 2
set key off
plot [50:200][:] \
"force.09.dat" u 2:($3/75) w l lw 2, \
"force.10.dat" u 2:($3/75) w l lw 2, \
"force.11.dat" u 2:($3/75) w l lw 2, \
"force.12.dat" u 2:($3/75) t "12", 9.7674460679523656e-01

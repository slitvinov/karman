set macro
u = 'u 2:($3/75) w l lw 2
set key off
plot [60:][0:] \
     "force.09.400.dat" @u, \
     "force.10.400.dat" @u, \
     "force.11.400.dat" @u, \
     "force.12.400.dat" @u

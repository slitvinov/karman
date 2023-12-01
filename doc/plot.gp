set term pngcairo
set output "nohelix.800.png"
set macro
u = 'u 2:($3/5) w l lw 2'
set key center top
set xlabel "time"
set ylabel "drag coeficient"

plot [0.2:][:2] \
"nohelix/force.08.800.dat" @u title "level: 08", \
"nohelix/force.09.800.dat" @u title "level: 09", \
"nohelix/force.10.800.dat" @u title "level: 10"

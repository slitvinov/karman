set term pngcairo
set output "nohelix.2000.png"
set macro
u = 'u 2:($3/5) w l lw 2'
set key center bottom
set xlabel "time"
set ylabel "drag coeficient"

plot [0.2:][0:] \
"nohelix/force.08.2000.dat" @u title "level: 08", \
"nohelix/force.09.2000.dat" @u title "level: 09", \
"nohelix/force.10.2000.dat" @u title "level: 10"

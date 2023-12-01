set term pngcairo
set output "smooth.800.png"
set macro
u = 'u 2:($3/5) w l lw 2'
set key center bottom
set xlabel "time"
set ylabel "drag coeficient"

plot [0.2:][0:] \
"smooth/force.08.800.dat" @u title "level: 08", \
"smooth/force.09.800.dat" @u title "level: 09", \
"smooth/force.10.800.dat" @u title "level: 10"

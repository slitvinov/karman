set term pngcairo size 1422, 597
set output "jiang.png"
set key center bottom
set xlabel "time"
set ylabel "C"
set macro
s = 'w l lw 2 lc rgb'
plot [:][-1:1.5] \
     0.9, \
     "<sort -k 2 -g 220/[0123].dat" u 2:($3/25) @s "#ff0000" t "", \
     "" u 2:($4/25) @s "#ff0000" t "Re = 220", \
     "<sort -k 2 -g 240/0.dat" u 2:($3/25) @s "#00ff00" t "", \
     "" u 2:($4/25) @s "#00ff00" t "Re = 240", \
     "<sort -k 2 -g 2000/0.dat" u 2:($3/25) @s "#0000ff" t "", \
     "" u 2:($4/25) @s "#0000ff" t "Re = 2000"


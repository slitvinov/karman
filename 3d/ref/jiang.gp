set term pngcairo size 1500, 500
set output "jiang.png"
# set key off
set key off
plot [:100][:] \
     "<sort -k 2 -g 220/[0123].dat" u 2:($3/25) w l lw 0.5 lc 0, "" u 2:($4/25) w l lw 0.5 lc 0 t "Re = 220", \
     "<sort -k 2 -g 240/0.dat" u 2:($3/25) w l lw 2 lc 1, "" u 2:($4/25) w l lw 2 lc 1 t "Re = 240", \
     "<sort -k 2 -g 2000/0.dat" u 2:($3/25) w l lw 2 lc 2, "" u 2:($4/25) w l lw 2 lc 2 t "Re = 2000"


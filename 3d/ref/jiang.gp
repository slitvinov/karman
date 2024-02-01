set term pngcairo size 1500, 500
set output "jiang.png"
set key off
plot [:4500][-1:1.5] "<cat [012].dat | sort -g" u 2:($3/25) w l lw 0.5, "" u 2:($4/25) w l lw 0.5

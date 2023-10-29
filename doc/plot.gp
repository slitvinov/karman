set macro
u = 'u 2:($3/75) w l lw 2
set key off
# set log y
plot [50:1400][0:] \
            "<sh post.sh" u 2:3:1 w p ps 2 pt 7, \
	    1.4

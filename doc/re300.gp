set macro
u = 'u 2:($3/75) w l lw 2
set key off
plot  [200:1200][0:] \
	    "<sh post.sh | awk '$1 == 10'" u 2:3, \
	    "../3d/ref/tritton.interp.txt"

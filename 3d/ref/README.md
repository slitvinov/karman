<pre>
for r in 2000 220 240
do for i in 0 1 2 3 4 5
   do rsync rc:/n/holyscratch01/koumoutsakos_lab/slitvinov/$r/$i/force.dat $r/$i.dat
   done
done
</pre>

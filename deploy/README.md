Compile

<pre>
$ make cylinder
</pre>
or build also MPI version
<pre>
$ make cylinder
</pre>

Help message

<pre>
$ ./cylinder -h
Usage: cylinder [-h] [-i] [-v] -r <Reynolds number> -l <resolution level> -p <dump period> -e <end time>
Options:
  -h     Display this help message
  -v     Verbose
  -i     Enable PPM image dumping
  -r <Reynolds number>     the Reynolds number (a decimal number)
  -l <resolution level>    the resolution level (positive integer)
  -p <dump period>         the dump period (positive integer)
  -e <end time>            end time of the simulation (decimal number)

Example usage:
  ./cylinder -v -i -r 100 -l 10 -p 100 -e 2
</pre>

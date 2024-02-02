Compile

<pre>
$ make cylinder
</pre>

Help message
<pre>
Usage: cylinder [-h] [-i] [-v] -r <Reynolds number> -l <resolution level> -p <dump period> -e <end time>
Options:
  -h     Display this help message
  -v     Verbose
  -i     Enable PPM image dumping
  -r <Reynolds number>     the Reynolds number (a decimal number)
  -l <resolution level>    the resolution level (positive integer)
  -o <preifx>              a prefix for the output files
  -p <dump period>         the dump period (positive integer)
  -e <end time>            end time of the simulation (decimal number)

Example usage:
  ./cylinder -v -i -r 100 -l 10 -p 100 -e 2 -o h
</pre>

<h2>Karman Vortex Street</h2>

This repository uses
[Basilisk](http://basilisk.fr/src/INSTALL)
for simulating Karman vortex street phenomena. Follow these steps for
installation:

<pre>
$ wget http://basilisk.fr/basilisk/basilisk.tar.gz
$ tar zxvf basilisk.tar.gz
$ cd basilisk/src
</pre>

On Linux:
<pre>
$ cp config.gcc config
</pre>

On macOS
<pre>
$ cp config.osx config
</pre>

<pre>
$ make qcc
$ cp qcc $HOME/.local/bin/
</pre>

See [deploy/README.md](deploy/README.md).

<h2>Results</h2>

<p align="center"><img src="img/karman.gif"/></p>

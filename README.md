# Karman Vortex Street

This repository uses
[Basilisk](http://basilisk.fr/src/INSTALL)
for simulating Karman vortex street phenomena. Follow these steps for
installation:

```
wget http://basilisk.fr/basilisk/basilisk.tar.gz
tar zxvf basilisk.tar.gz
cd basilisk/src
```

On Linux:
```
cp config.gcc config
```

On macOS
```
cp config.osx config
```

```
make
cp qcc $HOME/.local/bin/
```

See [deploy/README.md](deploy/README.md).

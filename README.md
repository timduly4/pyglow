![alt text](http://remote2.csl.illinois.edu/~duly/pyglow/logo.png "pyglow")

# Overview

pyglow is a Python module that wraps several upper atmosphere climatoglogical models written in FORTRAN.

# Installation

1. Download `./dist/pyglow-X.X.tar.gz` and unpack it in your local directory

2. Download the climatological models and wrap them with f2py:

```
$ cd ./pyglow/models/
$ make all
```

* If successful, there should be a `*.so` file in each of the `./models/dl_models/<model>/` directories

3. Install the python package

```
$ python ./setup.py install 
```

* On a mac, the folder `pyglow` and `*.so` files from `./models/dl_models/<model>/` should be in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`

* If you are denied permission, I recommend adding `--user` flag in command

4.  Test

* In Python, run:

```
from pyglow.pyglow import Point
from datetime import datetime

dn = datetime(2010, 3, 23, 9, 30)
lat = 0.
lon = -80.
alt = 250.

pt = Point(dn, lat, lon, alt)

pt.run_igrf()
pt.run_hwm93()
pt.run_msis()
pt.run_iri()

print pt
```

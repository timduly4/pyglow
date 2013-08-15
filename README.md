![alt text](https://raw.github.com/timduly4/pyglow/master/logo.png "pyglow")
[_(airglow viewed aboard the ISS)_](http://en.wikipedia.org/wiki/File:Cupola_above_the_darkened_Earth.jpg)

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
    $ cd ../../   # get back to root directory
    $ python ./setup.py install 
```
  * On a mac, the folder `pyglow` and `*.so` files from `./models/dl_models/<model>/` should be in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`
  * If you are denied permission, I recommend adding `--user` flag in command

# Testing / Example

* In Python, run:

```
from pyglow.pyglow import Point
from datetime import datetime

dn = datetime(2011, 3, 23, 9, 30)
lat = 0.
lon = -80.
alt = 250.

pt = Point(dn, lat, lon, alt)

print "Before running any models:"
print pt

pt.run_igrf()
pt.run_hwm93()
pt.run_msis()
pt.run_iri()

print "After running models:"
print pt
```

* The output should be as follows:

```
Before running any models:
               dn = 2011-03-23 09:30:00
  (lat, lon, alt) = (0.00, -80.00, 250.00)

Geophysical Indices:
---------------------
               kp = 2.70
               ap = 12.00
            f107a = 110.47

From IGRF:
-----------
     (Bx, By, Bz) = ( nan,  nan,  nan)
                B =  nan
              dip =  nan

From HWM:
------------
      hwm version = nan
           (u, v) = ( nan,  nan)

After running models:
               dn = 2011-03-23 09:30:00
  (lat, lon, alt) = (0.00, -80.00, 250.00)

Geophysical Indices:
---------------------
               kp = 2.70
               ap = 12.00
            f107a = 110.47

From IGRF:
-----------
     (Bx, By, Bz) = (2.448e-05, -6.449e-07, 9.932e-06)
                B = 2.642e-05
              dip = 22.08

From HWM:
------------
      hwm version = 93
           (u, v) = (14.07, 12.18)
```


# Hints

1. Use tab completion in ipython to view the full set of member data and variables available in the Point class.
  * For example, in the test code, run `pt.<TAB><TAB>` and class information will be listed.

# Uninstallation 

1. The install directory for pyglow is outputted when you run the `python ./setup.py install` command.  For example, on a mac this is usually in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`.
2.  Simply remove the `*.so` climatological models in this directory, as well as the `pyglow` and `pyglow_trash` folders.

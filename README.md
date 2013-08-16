![alt text](https://raw.github.com/timduly4/pyglow/master/logo.png "pyglow")
[_(airglow viewed aboard the ISS)_](http://en.wikipedia.org/wiki/File:Cupola_above_the_darkened_Earth.jpg)

# Overview

pyglow is a Python module that wraps several upper atmosphere climatoglogical models written in FORTRAN.

# Installation

1. Download `./dist/pyglow-X.X.tar.gz` and unpack it in your local directory
  * Alternatively, you can clone the repository:
```
    git clone git://github.com/timduly4/pyglow.git
```

2. Download the climatological models and wrap them with f2py:
```
    $ cd ./pyglow/models/
    $ make all
```
  * If successful, there should be a `*.so` file in each of the `./models/dl_models/<model>/` directories:

    ```
    $ find . -name "*.so"
    ./dl_models/hwm07/hwm07py.so
    ./dl_models/hwm93/hwm93py.so
    ./dl_models/igrf/igrf11py.so
    ./dl_models/iri/iri12py.so
    ./dl_models/msis/msis00py.so
    ```

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

## General
1. Use tab completion in ipython to view the full set of member data and variables available in the Point class.
  * For example, in the test code, run `pt.<TAB><TAB>` and class information will be listed.

## Updating geophysical indices
1. You'll need to download the geophysical indices as they become available.  The `update_kpap()` function is available in pyglow that enables you do this:

```
from pyglow.pyglow import update_kpap
update_kpap([2012, 2013]) # grabs indices for 2012 and 2013
update_kpap() # grabs all indices starting from 1932 to the current year
```
  * Only need to run this function when you would like to update the indices.
test


# Uninstallation 

1. The install directory for pyglow is outputted when you run the `python ./setup.py install` command.  For example, on a mac this is usually in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`.
2.  Simply remove the `*.so` climatological models in this directory, as well as the `pyglow` and `pyglow_trash` folders.

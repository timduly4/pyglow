![alt text](https://raw.github.com/timduly4/pyglow/master/logo.png "pyglow")
[_(airglow viewed aboard the ISS)_](http://en.wikipedia.org/wiki/File:Cupola_above_the_darkened_Earth.jpg)

# Overview

`pyglow` is a Python module that wraps several upper atmosphere climatoglogical models written in FORTRAN, such as the Horizontal Wind Model (HWM), the International Geomagnetic Reference Field (IGRF), the International Reference Ionosphere (IRI), and the Mass Spectrometer and Incoherent Scatter Radar (MSIS).

It includes the following upper atmospheric models:

  * HWM 2014
  * HWM 2007
  * HWM 1993
  * IGRF 12
  * IGRF 11
  * IRI 2016
  * IRI 2012
  * MSIS 2000

pyglow also provides access to the the following geophysical indices:
  * AP
  * Kp
  * F10.7
  * DST
  * AE

`pyglow` offers access to these models & indices in a convenient, high-level object-oriented interface within Python.

# Prerequisites

`pyglow` requires the following packages for installation:

1. `gfortran` (`$ sudo apt-get install gfortran`)
2. `ipdb` (`$ pip install ipdb --user`)

# Installation

### I'm Feeling Lucky:

First, checkout the repository:

```
    $ git clone git://github.com/timduly4/pyglow.git pyglow
```

Change directories into the repository folder and run the installation script:
```
    $ cd pyglow/
    $ ./pyglow_install.sh
```

### Individual installation steps:

If you have troubles, follow the individual installation steps:

(1) Download the package:
```
    $ git clone git://github.com/timduly4/pyglow.git
    $ cd pyglow/
```

(2) Download the climatological models and wrap them with f2py:
```
    $ cd ./pyglow/models/
    $ make all
```
  * If successful, there should be a `*.so` file in each of the `./models/dl_models/<model>/` directories:

    ```
    $ find . -name "*.so"
    ./dl_models/hwm07/hwm07py.so
    ./dl_models/hwm93/hwm93py.so
    ./dl_models/hwm14/hwm14py.so
    ./dl_models/igrf11/igrf11py.so
    ./dl_models/igrf12/igrf12py.so
    ./dl_models/iri12/iri12py.so
    ./dl_models/iri16/iri16py.so
    ./dl_models/msis/msis00py.so
    ```

(3) Install the python package
```
    $ cd ../../   # get back to root directory
    $ python ./setup.py install --user --prefix=
```
  * On a mac, the folder `pyglow` and `*.so` files from `./models/dl_models/<model>/` should be in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`
  * If you are denied permission, I recommend adding `--user` flag in command

(4) Download the geophysical indices

```
	$ cd ~/
	$ python -c "import pyglow; pyglow.update_indices()"
```


# Testing / Examples

See example scripts located in `./tests` for example calls to `pyglow`.

# Hints

### General
1. Use tab completion in ipython to view the full set of member data and variables available in the Point class.
  * For example, in the test code, run `pt.<TAB><TAB>` and class information will be listed.

### Updating geophysical indices with `update_indices()`
1. You'll need to download the geophysical indices as they become available.  The `update_indices()` function is available in pyglow that enables you do this:

```
from pyglow.pyglow import update_indices
update_indices([2012, 2013]) # grabs indices for 2012 and 2013
update_indices() # grabs all indices starting from 1932 to the current year
```

  * Only need to run this function when you would like to update the indices.


# Uninstallation 

1. The install directory for pyglow is outputted when you run the `python ./setup.py install` command.  For example, on a mac this is usually in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`.  You can also retrieve the library location with `python -c "import pyglow; print(pyglow.__file__)"`
2.  Simply remove the `*.so` climatological models in this directory, as well as the `pyglow` and `pyglow_trash` folders.


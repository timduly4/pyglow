		International Reference Ionosphere 2018 (IRI-2018) Software
		----------------------------------------------------------

================================================================================

The IRI is a joined project of the Committee on Space Research (COSPAR) and the 
International Union of Radio Science (URSI). IRI is an empirical model specifying
monthly averages of electron density, ion composition, electron temperature, and 
ion temperature in the altitude range from 50 km to 1500 km.

This directory includes the FORTRAN program, coefficients, and indices files for 
the latest version of the International Reference Ionosphere model: IRI2012. This
version includes several options for different parts and parameters. A logical
array JF(30) is used to set these options. The IRI-recommended set of options
can be found in the COMMENT section at the beginning of IRISUB.FOR. IRITEST.FOR 
sets these options as the default.

The compilation/link command in Fortran 77 is:
f77 -o iri iritest.for irisub.for irifun.for iritec.for iridreg.for igrf.for 
 cira.for iriflip.for

================================================================================
IMPORTANT: Please note that the following files have NOT been changed and are the
same as for IRI2007: CCIR%%, URSI%%
================================================================================
 
Directory Contents:
-----------------------------------------------------------------------------------

00_iri2018.tar	TAR file that includes all files from this directory and that was 
                created in UNIX using 'tar -cvf 00_iri2012.tar *'. UNIX command to 
                unpack is "tar -xvf 00_iri2012.tar". 
				 
irisub.for      This file contains the main subroutine iri_sub. It computes 
                height profiles if IRI output parameters (Ne, Te, Ti, Ni, vi) 
                for specified date and location. Also included is the 
                subroutine iri_web that computes output parameters versus
                latitude, longitude (geog. or geom.), year, day of year, hour 
                (LT or UT), and month. An example of how to call and use iri_web 
                is shown in iritest.for. Compilation of iritest.for requires
                irisub.for, irifun.for, iritec.for, iridreg.for, igrf.for, and
                cira.for.

irifun.for      This file contains the subroutines and functions that are 
                required for running IRI.


iriflip.for     Subroutines for the FLIP-related new model for the bottomside
                ion composition of Richards et al.
                
iridreg.for     Subroutines for the D region models of Friedrich-Torkar
                and of Danilov et al.

iritec.for      This file includes the subroutines for computing the ionospheric 
                electron content from 60km up to a specified upper limit. 

cira.for        This file includes the subroutines and functions for computing 
                the COSPAR International Reference Ionosphere (NRLMSIS-00) 
                neutral temperature and densities. 

igrf.for        This file includes the subroutines for the International
                Geomagnetic Reference Field (IGRF).

dgrf%%%%.dat    Definitive IGRF coefficients for the years 1945 to 2010 in steps
                of 5 years (%%%%=1945, 1950, etc.)(ASCII).
igrf%%%%.dat    Prelimenary IGRF coefficients for most recent year (ASCII).
igrf%%%%s.dat   IGRF coefficients for extrapolating 5 years into the future (ASCII).
 
MCSAT%%.dat		Monthly coefficient files for the Shubin(2015) COSMIC-based hmF2 model
                %%=month+10
                
iritest.for 	Test program indicating how to use of the IRI subroutines. 
                Requires irisub, irifun, iritec, iridreg, igrf,and cira.
                
ADDITIONAL FILES REQUIRED BUT NOT INCLUDED IN THIS DIRECTORY. THESE FILES ARE
AVAILABLE SEPARATELY ON THE IRI HOMEPAGE IRIMODEL.ORG.

ig_rz.dat       This file(s) contains the solar and ionospheric indices (IG12, Rz12) 
                for the time period from Jan 1958 onward. The file is updated 
                quarterly. It is read by subroutine tcon in irifun.for (ASCII). 
                [This file will be updated at close to quarterly intervals]
                
apf107.dat      This file provides the 3-hour ap magnetic index and F10.7 daily
                81-day and 365-day index from 1960 onward (ASCII).
                [This file will be updated at close to quarterly intervals]

CCIR%%.dat		Monthly coefficient files for the CCIR foF2 and M(3000)F2 models
                %%=month+10

URSI%%.dat		Monthly coefficient files for the URSI foF2 model
                %%=month+10

-----------------------------------------------------------------------------------
-----------------------------------------------------------------------------------

NOTE: Please consult the 'listing of changes' in the comments section at the top 
of each one of these programs for recent corrections and improvements.

More information about the IRI project can be found at  http://IRI.gsfc.nasa.gov

IRI parameters can be calculated and plotted online at sites given on the homepage


----------------------------- dbilitza@gmu.edu -------------------------------------

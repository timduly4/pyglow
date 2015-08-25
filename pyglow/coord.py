
from scipy import mat, cos, sin, arctan, sqrt, pi, arctan2, e
import numpy as np

""" GPS Constants"""
a = 6378137;                # semi-major axis of the earth [m]
b = 6356752.3145;           # semi-minor axis of the earth [m]
e = sqrt(1-(b**2)/(a**2));  # eccentricity of the earth = 0.08181919035596
lat_accuracy_thresh = 1.57e-6; #  10 meter latitude accuracy
leapSeconds = 14;           # GPS leap seconds [s]
muearth = 398600.5e9;       # G*Me, the "gravitational constant" for orbital
                            # motion about the Earth [m^3/s^2]
OmegaEDot = 7.2921151467e-5;# the sidereal rotation rate of the Earth
                            # (WGS-84) [rad/s]
c = 299792458;              # speed of light [m/s]
F = -4.442807633e-10;       # Relativistic correction term [s/m^(1/2)]
f1 = 1.57542e9;             # GPS L1 frequency [Hz]
f2 = 1.22760e9;             # GPS L2 frequency [Hz]

def say_hello():
    print("hello from coord.py!")

def ecef2lla(xyz):
    """
    Function: ecef2lla(xyz)
    ---------------------
    Converts ECEF X, Y, Z coordinates to WGS-84 latitude, longitude, altitude

    Inputs:
    -------
        xyz : 1x3 vector containing [X, Y, Z] coordinate
        

    Outputs:
    --------
        lla : 1x3 vector containing the converted [lat, lon, alt]
              (alt is in [m])  

    Notes:
    ------
        Based from Jonathan Makela's GPS_WGS84.m script

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
# TODO 
# [ ] make it vectorizable ?
    x = xyz[0][0]
    y = xyz[0][1]
    z = xyz[0][2]

    run = 1

    lla = np.array(np.zeros(xyz.size))
    # Compute longitude:
    lla[1] = arctan2(y,x)*(180./pi)

    # guess iniital latitude (assume you're on surface, h=0)
    p = sqrt(x**2+y**2)
    lat0 = arctan(z/p*(1-e**2)**-1)

    while (run == 1):
        # Use initial latitude to estimate N:
        N = a**2 / sqrt(a**2 * (cos(lat0))**2+b**2*(sin(lat0))**2)

        # Estimate altitude 
        h = p/cos(lat0)-N
    
        # Estimate new latitude using new height:
        lat1 = arctan(z/p*(1-((e**2*N)/(N+h)))**-1)
        
        if abs(lat1-lat0) < lat_accuracy_thresh:
            run = 0
        
        # Replace our guess latitude with most recent estimate:
        lat0 = lat1

    # load output array with best approximation of latitude (in degrees)
    # and altiude (in meters)
    
    lla[0] = lat1*(180./pi)
    lla[2] = h
    
    return lla


def lla2ecef(lla):
    """
    Function: lla2ecef(lla)
    ---------------------
    Converts WGS-84 latitude, longitude, altitude to ECEF X, Y, Z coordinates

    Inputs:
    -------
        lla : 1x3 vector containing the converted [lat, lon, alt]
              (alt is in [m])  

    Outputs:
    --------
        xyz : 1x3 vector containing [X, Y, Z] coordinate

    Notes:
    ------
        Based from Jonathan Makela's GPS_ECEF.m script

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)
        9/11/12 Updated to include vectorization.

    """
    # TODO 
    # [x] make it vectorizable ?
    
    # check for 1D case:
    dim = len(lla.shape)
    if dim == 1:
        lla = np.reshape(lla,(1,3))

    # convert lat and lon to radians
    lat = lla[:,0]/180.*pi
    lon = lla[:,1]/180.*pi
    alt = lla[:,2];

    xyz = np.array(np.zeros(lla.shape))

    N = a**2/sqrt((a*cos(lat))**2+(b*sin(lat))**2)

    # Calculate the X-coordinate
    xyz[:,0] = (N+alt)*cos(lat)*cos(lon)

    # Calculate the Y-coordinate
    xyz[:,1] = (N+alt)*sin(lon)*cos(lat)

    # Calculate the Z-coordinate
    xyz[:,2] = (N*(1-e**2)+alt)*sin(lat)

    return np.array(xyz)
    

def ven2ecef(lla,ven):
    """
    Function: ven2ecef(lla,ven)
    ---------------------
    Convert a vector given in VEN coordinates to ECEF coordinates

    Inputs:
    -------
        lla : 1x3 vector containing the converted [lat, lon, alt]
              (alt is in [m])  
        ven : 1x3 vector given [vertical, east, north] coordinates

    Outputs:
    --------
        xyz : 1x3 vector containing [X, Y, Z] coordinate

    Notes:
    ------
        Based from Jonathan Makela's GPS_VEN2ECEF.m script

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """

    # convert reference location to ECEF:
    ecef = lla2ecef(lla)


    #print ecef
    Xr = np.array(ecef[0][0])
    Yr = np.array(ecef[0][1])
    Zr = np.array(ecef[0][2])
       
    # convert to radians:
    refLong = pi/180.*lla[1]
  
    # calculate the geocentric latitude
    phiP = arctan2(Zr,sqrt(Xr**2+Yr**2))

    # calculate the ECEF location of the point    
    X = -sin(refLong)*ven[1] - cos(refLong)*sin(phiP)*ven[2]+cos(refLong)*cos(phiP)*ven[0]+Xr
    Y =  cos(refLong)*ven[1] - sin(refLong)*sin(phiP)*ven[2]+cos(phiP)*sin(refLong)*ven[0]+Yr
    Z = cos(phiP)*ven[2]+sin(phiP)*ven[0]+Zr

    # Subtract out the reference location
    XYZ = np.array([float(X-Xr), float(Y-Yr), float(Z-Zr)])
    
    return XYZ

def ecef2enu(XYZr, XYZp, lat_r, lon_r):
    return lat_r
    

def test_coord():
    lat = 30.
    lon = -10.
    alt = 300.*1e3

    print("ecef2lla( lla2ecef([%3.1f, %3.1f, %3.0e]) ) =" % (lat, lon, alt))
    print(ecef2lla( lla2ecef([lat, lon, alt]) ))
    
    xyz = np.array([-3197773.77194971,  -563853.79419661, -5587079.67459298])
    print("\nlla2ecef( ecef2lla([%3.6e, %3.6e, %3.6e]) ) =" % tuple(xyz))
    print(lla2ecef( ecef2lla(xyz) )) 

    lla = [0,0,300e3]
    ven = [0,0,200e3]
    print("\ntesting ven2ecef()...")
    print(ecef2lla(lla2ecef(lla) + ven2ecef(lla,ven)))
    


if __name__ == '__main__':
    test_coord()
#    [lat_r, lon_r, alt_r] = [20., 0., 280.e3]
#    [lat_p, lon_p, alt_p] = [20.5, 0., 270.e3]
#    XYZr = lla2ecef(np.array([lat_r, lon_r, alt_r]))
#    XYZp = lla2ecef(np.array([lat_p, lon_p, alt_p])) 

#    from numpy import sin, cos, pi
#    A = np.mat( [-sin(lat_r*pi/180.), cos(lat_r*pi/180.), 0.], \
#                [-sin(lon_r*pi/180.)*cos(lat_r*pi/180.), -sin(lon_r*pi/180.)*sin(lat_r*pi/180.), cos(lon_
    







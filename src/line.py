import numpy as np

from . import coord
from .point import Point


def _igrf_tracefield(dn, lat, lon, alt, target_ht, step):
    """
    Helper function to trace along a magnetic field line using IGRF

    :param dn: datetime.datetime object of requested trace
    :param lat: Latitude [degrees]
    :param lon: Longitude [degrees]
    :param alt: Altitude [km]
    :param target_ht: Altitude to stop trace [km]
    :param step: Step size of trace [km]

    :return lla: (latitude, longitude, altitude) data structure of trace
                 with dimentions [Nsteps x 3]

    """

    # Go North:
    lla_north = _igrf_tracefield_hemis(
        dn,
        lat,
        lon,
        alt,
        target_ht,
        step,
    )

    # Go South:
    lla_south = _igrf_tracefield_hemis(
        dn,
        lat,
        lon,
        alt,
        target_ht,
        -step,
    )

    # Stack them together:
    lla = np.vstack(
        [
            np.flipud(lla_north)[:-1, :],
            lla_south,
        ]
    )

    return lla


def _igrf_tracefield_hemis(dn, lat, lon, alt, target_ht, step):
    """
    Helper function to trace along a magnetic field line using IGRF
    for only one hemisphere

    :param dn: datetime.datetime object of requested trace
    :param lat: Latitude [degrees]
    :param lon: Longitude [degrees]
    :param alt: Altitude [km]
    :param target_ht: Altitude to stop trace [km]
    :param step: Step size of trace [km]

    :return lla: (latitude, longitude, altitude) data structure of trace
                 with dimentions [Nsteps x 3]

    """

    lat = float(lat)
    lon = float(lon)
    alt = float(alt)
    target_ht = float(target_ht)
    step = float(step)

    target_ht = target_ht*1e3
    step = step*1e3

    lla = np.array([lat, lon, alt*1e3])

    lla_field = lla

    """ Step 1: trace the field along a given direction """
    TOLERANCE = 10  # [m]
    i = 0
    while (lla[2] > target_ht):
        # convert to ECEF:
        ecef = coord.lla2ecef(lla)

        # Grab field line information:
        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        # coordinates follow pyglow's convention:
        # (x -> east, y -> north, z -> up)

        N = p.By  # North
        E = p.Bx  # East
        D = -p.Bz  # Down
        A = p.B   # Total

        # Step along the field line
        ecef_new = ecef + coord.ven2ecef(
            lla,
            [-D/A*step, E/A*step, N/A*step],
        )

        # Convert to lla coordinates:
        lla = coord.ecef2lla(ecef_new)

        # add the field line to our collection:
        lla_field = np.vstack([lla_field, lla])
        i = i + 1

    """ Step 2: Make the last point close to target_ht """
    while (abs(lla[2]-target_ht) > TOLERANCE):

        # Find out how much we need to step by:
        step = -np.sign(step)*abs(lla[2]-target_ht)

        ecef = coord.lla2ecef(lla)

        p = Point(dn, lla[0], lla[1], lla[2]/1e3)
        p.run_igrf()

        N = p.Bx  # North
        E = p.By  # East
        D = p.Bz  # Down
        A = p.B   # Total

        # Trace the field, but use the modified step:
        ecef_new = ecef + coord.ven2ecef(
            lla,
            np.array([-D/A, E/A, N/A]) * step/(-D/A),
        )

        # TODO : I changed this, is this correct?
        lla = coord.ecef2lla(ecef_new)

    # replace last entry with the point close to target_ht:
    lla_field[-1, :] = lla

    return lla_field


def Line(dn, lat, lon, alt, target_ht=90., step=15.):
    """
    Return a list of instances of Point by
    tracing along the geomagnetic field line.

    pts = Line(dn, lat, lon, alt, target_ht=90., step=15.)

    :param dn: datetime.datetime object of requested trace
    :param lat: Latitude [degrees]
    :param lon: Longitude [degrees]
    :param alt: Altitude [km]
    :param target_ht: (optional) Altitude to stop trace [km]
    :param step: (optional) Step size of trace [km]

    """
    llas = _igrf_tracefield(dn, lat, lon, alt, target_ht, step)
    pts = []

    for lla in llas:
        pts.append(
            Point(dn, lla[0], lla[1], lla[2]/1e3)
        )

    return pts

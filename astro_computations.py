import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.visualization import quantity_support


def body_altitude_curve(target_cord, telescope_location, date, delta_midnight=None):

    # target_name = 'M33'
    # telescope_location = dict(lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m)
    # date = "2012-7-13 00:00:00"

    # target_cord = SkyCoord.from_name(target_name)
    # m33 = SkyCoord(23.46206906, 30.66017511, unit="deg")
    # bear_mountain = EarthLocation(**telescope_location)

    utcoffset = -4 * u.hour  # EDT
    midnight = Time(date) - utcoffset

    times_July12_to_13 = midnight + delta_midnight
    frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=telescope_location)

    sunaltazs_trajectory = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
    moonaltazs_trajectory = get_body("moon", times_July12_to_13).transform_to(frame_July12_to_13)
    targetaltazs_trajectory = target_cord.transform_to(frame_July12_to_13)

    return sunaltazs_trajectory, moonaltazs_trajectory, targetaltazs_trajectory



# def body_altitude_curve(target_name=None, telescope_location=None, date=None, delta_midnight=None):
#
#     target_name = 'M33'
#     telescope_location = dict(lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m)
#     date = "2012-7-13 00:00:00"
#
#     target_cord = SkyCoord.from_name(target_name)
#     # m33 = SkyCoord(23.46206906, 30.66017511, unit="deg")
#     bear_mountain = EarthLocation(**telescope_location)
#
#     utcoffset = -4 * u.hour  # EDT
#     midnight = Time(date) - utcoffset
#
#     times_July12_to_13 = midnight + delta_midnight
#     frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bear_mountain)
#
#     sunaltazs_trajectory = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
#     moonaltazs_trajectory = get_body("moon", times_July12_to_13).transform_to(frame_July12_to_13)
#     targetaltazs_trajectory = target_cord.transform_to(frame_July12_to_13)
#
#     return sunaltazs_trajectory, moonaltazs_trajectory, targetaltazs_trajectory


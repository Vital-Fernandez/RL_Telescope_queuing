import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body, get_sun
from astropy.time import Time
from astropy.visualization import quantity_support

# m33 = SkyCoord.from_name("M33")
m33 = SkyCoord(23.46206906, 30.66017511, unit="deg")
bear_mountain = EarthLocation(lat=41.3 * u.deg, lon=-74 * u.deg, height=390 * u.m)

utcoffset = -4 * u.hour  # EDT
time = Time("2012-7-12 23:00:00") - utcoffset

m33altaz = m33.transform_to(AltAz(obstime=time, location=bear_mountain))
print(f"M33's Altitude = {m33altaz.alt:.2}")

midnight = Time("2012-7-13 00:00:00") - utcoffset

delta_midnight = np.linspace(-12, 12, 1000) * u.hour
times_July12_to_13 = midnight + delta_midnight
frame_July12_to_13 = AltAz(obstime=times_July12_to_13, location=bear_mountain)

sunaltazs_July12_to_13 = get_sun(times_July12_to_13).transform_to(frame_July12_to_13)
moonaltazs_July12_to_13 = get_body("moon", times_July12_to_13).transform_to(frame_July12_to_13)
m33altazs_July12_to_13 = m33.transform_to(frame_July12_to_13)

with quantity_support():

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    ax.plot(delta_midnight, sunaltazs_July12_to_13.alt, color="r", label="Sun")
    ax.plot(delta_midnight, moonaltazs_July12_to_13.alt, color=[0.75] * 3, ls="--", label="Moon")

    mappable = ax.scatter(delta_midnight, m33altazs_July12_to_13.alt, c=m33altazs_July12_to_13.az.value,
                          label="M33", lw=0, s=8, cmap="viridis")

    ax.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg, sunaltazs_July12_to_13.alt < (-0 * u.deg),
                    color="0.5", zorder=0)

    ax.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg, sunaltazs_July12_to_13.alt < (-18 * u.deg),
                    color="k", zorder=0)

    fig.colorbar(mappable).set_label("Azimuth [deg]")
    ax.legend(loc="upper left")
    ax.set_xlim(-12 * u.hour, 12 * u.hour)
    ax.set_xticks((np.arange(13) * 2 - 12) * u.hour)
    ax.set_ylim(0 * u.deg, 90 * u.deg)
    ax.set_xlabel("Hours from EDT Midnight")
    ax.set_ylabel("Altitude [deg]")
    ax.grid(visible=True)
    plt.show()
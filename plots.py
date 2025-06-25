import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import quantity_support
from astropy import units as u


def altitude_plotter(delta_midnight, sun_trajectory, moon_trajectory, target_trajectory):

    with quantity_support():

        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax.plot(delta_midnight, sun_trajectory.alt, color="r", label="Sun")
        ax.plot(delta_midnight, moon_trajectory.alt, color=[0.75] * 3, ls="--", label="Moon")

        mappable = ax.scatter(delta_midnight, target_trajectory.alt, c=target_trajectory.az.value,
                              label="M33", lw=0, s=8, cmap="viridis")

        ax.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg, sun_trajectory.alt < (-0 * u.deg),
                        color="0.5", zorder=0)

        ax.fill_between(delta_midnight, 0 * u.deg, 90 * u.deg, sun_trajectory.alt < (-18 * u.deg),
                        color="k", zorder=0)

        fig.colorbar(mappable).set_label("Azimuth [deg]")
        ax.legend(loc="upper left")
        ax.set_xlim(-12 * u.hour, 12 * u.hour)
        ax.set_xticks((np.arange(13) * 2 - 12) * u.hour)
        ax.set_ylim(0 * u.deg, 90 * u.deg)
        ax.set_xlabel("Hours from EDT Midnight")
        ax.set_ylabel("Altitude [deg]")
        ax.grid(visible=True)
        # plt.show()

    return fig
import numpy as np
import streamlit as st

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation

from plots import altitude_plotter
from astro_computations import body_altitude_curve

list_telescopes = EarthLocation.get_site_names()

st.title("Reinforcement learning for pptimized astronomical target scheduling")

st.write("This site allows users to input a list of astronomical targets, select a ground-based observatory, "
         "and specify a target observation date. It uses reinforcement learning to dynamically optimize target"
         " selection, aiming to maximize altitude and overall observability under varying sky conditions. "
         "The system learns to prioritize targets that offer the best scientific return within the constraints of"
         " visibility and telescope scheduling.")

with st.form("observation_form"):
    col1, col2 = st.columns(2)

    with col1:
        target_input = st.text_input("Target name or coordinates", placeholder="e.g. M31 or 10h21m00s +20d00m00s")

    with col2:
        telescope_name = st.selectbox("Telescope location", list(list_telescopes))
        obs_date = st.date_input("Observation date")

    submitted = st.form_submit_button("Submit")

if submitted:

    skycoord_obj = None
    try:
        skycoord_obj = SkyCoord.from_name(target_input)
        st.success(f"Resolved name to coordinates: RA={skycoord_obj.ra.deg:.3f}, Dec={skycoord_obj.dec.deg:.3f}")
    except Exception:
        try:
            skycoord_obj = SkyCoord(target_input, unit=(u.hourangle, u.deg))
            st.success(f"Parsed coordinates: RA={skycoord_obj.ra.deg:.3f}, Dec={skycoord_obj.dec.deg:.3f}")
        except Exception as e:
            st.error(f"Could not resolve input: {e}")

    telescope_location = EarthLocation.of_site(telescope_name)
    formatted_time = f"{obs_date.year}-{obs_date.month}-{obs_date.day} 00:00:00"
    obs_time = Time(formatted_time)

    if skycoord_obj:

        st.write('***')
        st.subheader("Input observations summary")

        geo_str = f"Lat: {telescope_location.to_geodetic().lat:.4f}, Lon: {telescope_location.to_geodetic().lon:.4f}, Height: {telescope_location.to_geodetic().height.to(u.m):.0f} m"
        tar_str = skycoord_obj.to_string('hmsdms', precision=2)

        st.write("📍 **SkyCoord:**", tar_str )
        st.write("🌍 **EarthLocation:**", geo_str)
        st.write("🕛 **Observation Time:**", formatted_time)

        st.write('***')
        st.subheader("Object visibility curves")

        delta_midnight = np.linspace(-12, 12, 1000) * u.hour
        sun_trajectory, moon_trajectory, target_trajectory = body_altitude_curve(skycoord_obj,
                                                                                 telescope_location,
                                                                                 obs_time,
                                                                                 delta_midnight=delta_midnight)
        fig_trajectory = altitude_plotter(delta_midnight, sun_trajectory, moon_trajectory, target_trajectory)

        st.pyplot(fig_trajectory)

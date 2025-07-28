# Description: This script is used to analyze the orbit model of an NGSO constellation.
# The script calculates the number of visible satellites from a ground station and the elevation angles of the satellites.
# The script uses the OrbitModel class from the sharc.satellite.ngso.orbit_model module to calculate the satellite positions.
import numpy as np
from pathlib import Path
import plotly.graph_objects as go

from sharc.parameters.parameters_mss_d2d import ParametersMssD2d
from sharc.satellite.ngso.orbit_model import OrbitModel
from sharc.satellite.utils.sat_utils import calc_elevation


if __name__ == "__main__":

    # Input parameters
    local_dir = Path(__file__).parent.resolve()
    param_file = local_dir / ".." / "input" / "parameters_mss_d2d_to_imt_dl_co_channel_system_A.yaml"
    params = ParametersMssD2d()
    params.load_parameters_from_file(param_file)
    orbit_params = params.orbits[0]
    print("Orbit parameters:")
    print(orbit_params)

    # Ground station location
    GROUND_STA_LAT = -15.7801
    GROUND_STA_LON = -42.9292
    MIN_ELEV_ANGLE_DEG = 5.0  # minimum elevation angle for visibility

    # Time duration in days for the linear time simulation
    TIME_DURATION_HOURS = 72

    # Random samples
    N_DROPS = 50000

    # Random seed for reproducibility
    SEED = 6

    # Instantiate the OrbitModel
    orbit = OrbitModel(
        Nsp=orbit_params.sats_per_plane,
        Np=orbit_params.n_planes,
        phasing=orbit_params.phasing_deg,
        long_asc=orbit_params.long_asc_deg,
        omega=orbit_params.omega_deg,
        delta=orbit_params.inclination_deg,
        hp=orbit_params.perigee_alt_km,
        ha=orbit_params.apogee_alt_km,
        Mo=orbit_params.initial_mean_anomaly
    )

    # Show visible satellites from ground-station
    total_periods = int(TIME_DURATION_HOURS * 3600 / orbit.orbital_period_sec)
    print(f"Total periods: {total_periods}")
    pos_vec = orbit.get_satellite_positions_time_interval(initial_time_secs=0,
                                                          interval_secs=10,
                                                          n_periods=total_periods)
    sat_altitude_km = orbit.apogee_alt_km  # altitude of the satellites in kilometers
    num_of_visible_sats_per_drop = []
    elev_angles = calc_elevation(GROUND_STA_LAT, pos_vec['lat'], GROUND_STA_LON, pos_vec['lon'], sat_altitude_km)
    elevation_angles_per_drop = elev_angles[np.where(np.array(elev_angles) > MIN_ELEV_ANGLE_DEG)]
    num_of_visible_sats_per_drop = np.sum(np.array(elev_angles) > MIN_ELEV_ANGLE_DEG, axis=0)

    # Show visible satellites from ground-station - random
    num_of_visible_sats_per_drop_rand = []
    elevation_angles_per_drop_rand = []
    rng = np.random.RandomState(seed=SEED)
    pos_vec = orbit.get_orbit_positions_random(rng=rng, n_samples=N_DROPS)
    elev_angles = calc_elevation(GROUND_STA_LAT,
                                 pos_vec['lat'],
                                 GROUND_STA_LON, pos_vec['lon'],
                                 sat_altitude_km)
    elevation_angles_per_drop_rand = elev_angles[np.where(np.array(elev_angles) > MIN_ELEV_ANGLE_DEG)]
    num_of_visible_sats_per_drop_rand = np.sum(np.array(elev_angles) > MIN_ELEV_ANGLE_DEG, axis=0)

    # Free some memory
    del elev_angles
    del pos_vec

    # plot histogram of visible satellites
    fig = go.Figure(data=[go.Histogram(x=num_of_visible_sats_per_drop,
                                       histnorm='probability', xbins=dict(start=-0.5, size=1))])
    fig.update_layout(
        title_text='Visible satellites per drop',
        xaxis_title_text='Num of visible satellites',
        yaxis_title_text='Probability',
        bargap=0.2,
        bargroupgap=0.1,
        xaxis=dict(
            tickmode='linear',
            tick0=0,
            dtick=1
        )
    )

    fig.add_trace(go.Histogram(x=num_of_visible_sats_per_drop_rand,
                               histnorm='probability',
                               xbins=dict(start=-0.5, size=1)))

    fig.data[0].name = 'Linear time'
    fig.data[1].name = 'Random'
    fig.update_layout(legend_title_text='Observation Type')
    file_name = Path(__file__).parent / "visible_sats_per_drop.html"
    fig.write_html(file=file_name, include_plotlyjs="cdn", auto_open=False)

    # plot histogram of elevation angles
    fig = go.Figure(data=[go.Histogram(x=np.array(elevation_angles_per_drop).flatten(),
                                       histnorm='probability', xbins=dict(start=0, size=5))])
    fig.update_layout(
        title_text='Elevation angles',
        xaxis_title_text='Elevation angle [deg]',
        yaxis_title_text='Probability',
        bargap=0.2,
        bargroupgap=0.1
    )
    file_name = Path(__file__).parent / "elevation_angles.html"
    fig.write_html(file=file_name, include_plotlyjs="cdn", auto_open=False)

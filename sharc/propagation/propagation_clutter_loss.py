# -*- coding: utf-8 -*-

from sharc.propagation.propagation import Propagation
from sharc.support.enumerations import StationType

import numpy as np
import scipy
import matplotlib.pyplot as plt


class PropagationClutterLoss(Propagation):

    """
        This Recommendation (ITU_R_P_2108-1/2021) provides methods for
        estimating loss through clutter at frequencies between 30 MHz and 100 GHz.

        The ITU Radiocommunication Assembly,
        considering
             a) that, for system planning and interference assessment it may be
             necessary to account for the attenuation suffered by radio waves in
             passing over or between buildings;
             b) that, where a terrestrial station may be shielded by buildings a
             detailed calculation for a general case can be difficult to formulate
             and losses due to clutter must be considered dependant on the
             deployment scenario;
             c) that, where terrestrial stations are in motion the clutter envir-
             onment of the radio path will be variable,
        recognizing
             a) that Recommendation ITU-R P.1411 contains data and models for
             short-range radio system, mainly within an urban environment from
             300 MHz to 100 GHz;
             b) that Recommendation ITU-R P.2040 contains basic expressions for
             reflection from and penetration through building materials, and a
             harmonised representation of building material electrical properties
             above about 100 MHz;
             c) that Recommendation ITU-R P.452 contains a prediction method for
             the evaluation of interference between stations on the surface of
             the Earth at frequencies from about 0.1 GHz to 50 GHz, accounting
             for both clear-air and hydrometeor scattering interference mechanisms;
             d) that Recommendation ITU-R P.1812 describes a propagation predict-
             ion method suitable for terrestrial point-to-area services in the
             frequency range 30 MHz to 6 000 MHz;
             e) that Recommendation ITU-R P.833 presents several models to enable
             the user to evaluate the effect of vegetation on radiowave signals
             between 30 MHz and 60 GHz;
             f) that Recommendation ITU-R P.2109 provides a statistical model
             for building entry loss for frequencies between about 80 MHz and 100 GHz,
        recommends
             that the material in Recomendation ITU-R P.2108-1/2021 be used to
             estimate clutter loss.
    """

    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates clutter loss according to Recommendation P.2108-0

        Parameters
        ----------
            distance (np.array) : distances between stations [m]
            frequency (np.array) : center frequency [MHz]
            elevation (np.array) : elevation angles [deg]
            loc_percentage (np.array) : Percentage locations range [0, 1[
                                        "RANDOM" for random percentage (Default = RANDOM)
            station_type (StationType) : if type is IMT or FSS_ES, assume terrestrial
                terminal within the clutter (ref § 3.2); otherwise, assume that
                one terminal is within the clutter and the other is a satellite,
                aeroplane or other platform above the surface of the Earth.

        Returns
        -------
            array with clutter loss values with dimensions of distance

        """
        f = kwargs["frequency"]
        loc_per = kwargs.pop("loc_percentage", "RANDOM")
        type = kwargs["station_type"]
        d = kwargs["distance"]

        if f.size == 1:
            f = f * np.ones(d.shape)

        if isinstance(loc_per, str) and loc_per.upper() == "RANDOM":
            p1 = self.random_number_gen.random_sample(d.shape)
            p2 = self.random_number_gen.random_sample(d.shape)
        else:
            p1 = loc_per * np.ones(d.shape)
            p2 = loc_per * np.ones(d.shape)

        if type is StationType.IMT_BS or type is StationType.IMT_UE or type is StationType.FSS_ES:
            clutter_type = kwargs["clutter_type"]
            if clutter_type == 'one_end':
                loss = self.get_terrestrial_clutter_loss(f, d, p1, True)
            else:
                loss = self.get_terrestrial_clutter_loss(f, d, p1, True) + self.get_terrestrial_clutter_loss(f, d, p2, False)
        else:
            theta = kwargs["elevation"]
            earth_station_height = kwargs["earth_station_height"]
            mean_clutter_height = kwargs["mean_clutter_height"]
            below_rooftop = kwargs["below_rooftop"]
            loss = self.get_spacial_clutter_loss(f, theta, p1, earth_station_height, mean_clutter_height)
            mult_1 = np.zeros(d.shape)
            num_ones = int(np.round(mult_1.size * below_rooftop / 100))
            indices = np.random.choice(mult_1.size, size=num_ones, replace=False)
            mult_1.flat[indices] = 1
            loss *= mult_1
        return loss

    def get_spacial_clutter_loss(
        self, frequency: float,
        elevation_angle: float,
        loc_percentage,
        earth_station_height,
        mean_clutter_height
    ):
        """
        Computes clutter loss according to ITU-R P.2108-2 §3.3 for Earth-space and aeronautical paths.

        Parameters:
        - frequency: Frequency (GHz), 0.5 <= frequency <= 100
        - elevation_angle: Elevation angle (degrees), 0 <= elevation_angle <= 90
        - loc_percentage: Percentage of locations (%), 0 < loc_percentage < 100
        - earth_station_height: Ground station height (m), >= 1
        - mean_clutter_height: Median clutter height (m): Low-rise: <= 8; Mid-rise: 8 < ... <= 20; High-rise: > 20

        Returns:
        - Lces: Clutter loss according to P.2108 §3.3 (dB)
        """
        # Converting to GHz
        frequency = frequency / 1000
        # --- Table 7: plos parameters ---
        if mean_clutter_height == "Low":
            ak, bk, ck = 4.9, 6.7, 2.6
            aC, bC = 0.19, 0
            aV, bV = 1.4, 74
        elif mean_clutter_height == "Mid":
            ak, bk, ck = -2.6, 6.6, 2.0
            aC, bC = 0.42, -6.7
            aV, bV = 0.15, 97
        else:
            ak, bk, ck = 2.4, 7.0, 1.0
            aC, bC = 0.19, -2.7
            aV, bV = 0.15, 98

        # --- Table 8: p_fclo_los parameters ---
        if mean_clutter_height == "Low":
            akp, bkp = 6.0, 0.07
            aCp, bC1p, bC2p = 0.15, 5.4, -0.3
            bC3p, bC4p = 3.2, 0.07
            cCp, aVp, bVp = -27, 1.6, -17
        elif mean_clutter_height == "Mid":
            akp, bkp = 3.6, 0.05
            aCp, bC1p, bC2p = 0.17, 13, -0.2
            bC3p, bC4p = 3.7, 0.05
            cCp, aVp, bVp = -41, 1, -21
        else:
            akp, bkp = 5.0, 0.003
            aCp, bC1p, bC2p = 0.17, 32.6, 0.012
            bC3p, bC4p = -23.9, -0.07
            cCp, aVp, bVp = -41, 1, -18

        # --- LoS probability (eqs. 7-8) ---
        Vmax = np.minimum(aV * earth_station_height + bV, 100)
        Ce = aC * earth_station_height + bC
        k = ((earth_station_height + ak) / bk) ** ck
        num = 1 - np.exp(-k * (elevation_angle + Ce) / 90)
        den = 1 - np.exp(-k * (90 + Ce) / 90)
        pLoS = np.maximum(0, Vmax * (num / den))

        # --- Conditional probability of Fresnel clearance (eqs. 9-10) ---
        Vmaxp = np.minimum(aVp * earth_station_height + bVp, 0) * frequency ** (-0.55) + 100
        Cep = (frequency * 1e9) ** aCp + bC1p * np.exp(bC2p * earth_station_height) + bC3p * np.exp(bC4p *
                                                                                                    earth_station_height) + cCp
        kp = akp * np.exp(bkp * earth_station_height)
        nump = 1 - np.exp(-kp * (elevation_angle + Cep) / 90)
        denp = 1 - np.exp(-kp * (90 + Cep) / 90)
        pFcLoS_LoS = np.maximum(0, Vmaxp * (nump / denp))
        pFcLoS = pLoS * pFcLoS_LoS / 100

        # Table 9 constants
        alpha1, alpha2 = 8.54, 0.056
        beta1, beta2 = 17.57, 6.32
        gamma1, gamma2 = 0.63, 0.19

        mu = alpha1 + beta1 * np.log(1 + (90 - elevation_angle) / 90) + frequency ** gamma1
        sigma = alpha2 + beta2 * np.log(1 + (90 - elevation_angle) / 90) + frequency ** gamma2

        # Table 10 constants
        theta1, theta2 = -3.542, -155.1
        sigma1, sigma2 = -1.06, -6.342

        # Prepare output array
        Lces = np.zeros_like(frequency, dtype=float)

        # Branch logic: mask arrays for each regime
        mask_nlos = loc_percentage > pLoS
        mask_fclo = loc_percentage <= pFcLoS
        mask_between = (~mask_nlos) & (~mask_fclo)

        # NLoS branch
        if np.any(mask_nlos):
            pp = (loc_percentage[mask_nlos] - pLoS[mask_nlos]) / (100 - pLoS[mask_nlos])
            Finv = -np.sqrt(2) * scipy.special.erfcinv(2 * pp)
            Lces_nlos = mu[mask_nlos] + sigma[mask_nlos] * Finv
            Lces[mask_nlos] = np.maximum(Lces_nlos, 6)

        # Fresnel-clear branch
        if np.any(mask_fclo):
            ratio = loc_percentage[mask_fclo] / pFcLoS[mask_fclo]
            Lces[mask_fclo] = (theta1 * np.exp(theta2 * ratio) + sigma1 * np.exp(sigma2 * ratio))

        # Between
        if np.any(mask_between):
            Lces[mask_between] = (
                6 * (loc_percentage[mask_between] - pFcLoS[mask_between]) /
                (pLoS[mask_between] - pFcLoS[mask_between])
            )

        return Lces

    def get_terrestrial_clutter_loss(
        self,
        frequency: float,
        distance: float,
        loc_percentage: float,
        apply_both_ends=True,
    ):
        """
        This method gives a statistical distribution of clutter loss. The model
        can be applied for urban and suburban clutter loss modelling. An
        additional "loss" is calculated which can be added to the transmission
        loss or basic transmission loss. Clutter loss will vary depending on
        clutter type, location within the clutter and movement in the clutter.
        If the transmission loss or basic transmission loss has been calculated
        using a model (e.g. Recommendation ITU-R P.1411) that inherently
        accounts for clutter over the entire path then the method below should
        not be applied. The clutter loss not exceeded for loc_percentage% of
        locations for the terrestrial to terrestrial path, "loss"", is given
        by this method. The clutter loss must not exceed a maximum value
        calculated for 𝑑 = 2 𝑘m (loss_2km)
.

        Parameters
        ----
            frequency : center frequency [MHz] - Frequency range: 0.5 to 67 GHz
            distance : distance [m] - Minimum path length: 0.25 km (for the
                       correction to be applied at only one end of the path)
                                                           1.0 km (for the
                       correction to be applied at both ends of the path)
            loc_percentage : percentage of locations [0,1] - Percentage
                             locations range: 0 < p < 100

            apply_both_ends : if correction will be applied at both ends of the
                              path

        Returns
        -------
            loss : The clutter loss not exceeded for loc_percentage% of
                   locations for the terrestrial to terrestrial path
        """

        d = distance.copy()
        d = d.reshape((-1, 1))
        f = frequency.reshape((-1, 1))
        p = loc_percentage.reshape((-1, 1))

        sigma_l = 4.0
        sigma_s = 6.0

        loss = np.zeros(d.shape)
        loss_2km = np.zeros(d.shape)

        if apply_both_ends:
            # minimum path length for the correction to be applied at only one end of the path
            id_d = np.where(d >= 1000)[0]
        else:
            # minimum path length for the correction to be applied at both ends of the path
            id_d = np.where(d >= 250)[0]

        if len(id_d):
            Ll = -2.0 * \
                np.log10(
                    10 ** (-5.0 * np.log10(f[id_d] * 1e-3) - 12.5) + 10 ** (-16.5))
            Ls_temp = 32.98 + 3.0 * np.log10(f[id_d] * 1e-3)
            Ls = 23.9 * np.log10(d[id_d] * 1e-3) + Ls_temp
            invQ = np.sqrt(2) * scipy.special.erfcinv(2 * (p[id_d]))
            sigma_cb = np.sqrt(((sigma_l**(2.0)) * (10.0**(-0.2 * Ll)) + (sigma_s**(2.0)) *
                               (10.0**(-0.2 * Ls))) / (10.0**(-0.2 * Ll) + 10.0**(-0.2 * Ls)))
            loss[id_d] = -5.0 * \
                np.log10(10 ** (-0.2 * Ll) + 10 **
                         (-0.2 * Ls)) - sigma_cb * invQ

            # The clutter loss must not exceed a maximum value calculated for 𝑑 = 2 𝑘m (loss_2km)
            Ls_2km = 23.9 * np.log10(2) + Ls_temp
            sigma_cb_2km = np.sqrt(((sigma_l**(2.0)) * (10.0**(-0.2 * Ll)) + (sigma_s**(2.0)) *
                                   (10.0**(-0.2 * Ls_2km))) / (10.0**(-0.2 * Ll) + 10.0**(-0.2 * Ls_2km)))
            loss_2km[id_d] = -5.0 * \
                np.log10(10 ** (-0.2 * Ll) + 10 ** (-0.2 * Ls_2km)) - \
                sigma_cb_2km * invQ
            id_max = np.where(loss >= loss_2km)[0]
            loss[id_max] = loss_2km[id_max]

        loss = loss.reshape(distance.shape)

        return loss


if __name__ == '__main__':

    elevation_angle = np.array([90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5, 0])
    loc_percentage = np.linspace(0.1, 99.9, 1001)
    frequency = 3000 * np.ones(loc_percentage.shape)
    earth_station_height = 5 * np.ones(loc_percentage.shape)
    random_number_gen = np.random.RandomState(101)
    cl = PropagationClutterLoss(random_number_gen)
    clutter_loss = np.empty([len(elevation_angle), len(loc_percentage)])

    for i in range(len(elevation_angle)):
        clutter_loss[i, :] = cl.get_spacial_clutter_loss(
            frequency,
            elevation_angle[i],
            loc_percentage,
            earth_station_height,
            'Low',
        )

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')
    ax = fig.gca()

    for j in range(len(elevation_angle)):
        ax.plot(clutter_loss[j, :], loc_percentage,
                label="%i deg" % elevation_angle[j], linewidth=1)

    plt.title("Cumulative distribution of clutter loss not exceeded for 30 GHz")
    plt.xlabel("clutter loss [dB]")

    plt.ylabel("percent of locations [%]")
    plt.xlim((-5, 30))
    plt.ylim((0, 100))
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.grid()
    plt.show()

    distance = np.linspace(250, 100000, 100000)
    frequency = np.array([1, 2, 4, 8, 16, 32, 67]) * 1e3

    loc_percentage = 0.5 * np.ones(distance.shape)
    apply_both_ends = False

    clutter_loss_ter = np.empty([len(frequency), len(distance)])

    for i in range(len(frequency)):
        clutter_loss_ter[i, :] = cl.get_terrestrial_clutter_loss(
            frequency[i] * np.ones(distance.shape),
            distance,
            loc_percentage,
            apply_both_ends,
        )

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    # ax.set_prop_cycle( cycler('color', ['k', 'r', 'b', 'g']) )

    for j in range(len(frequency)):
        freq = frequency[j] * 1e-3
        ax.semilogx(distance * 1e-3,
                    clutter_loss_ter[j, :], label="%i GHz" % freq, linewidth=1)

    plt.title("Median clutter loss for terrestrial paths")
    plt.xlabel("Distance [km]")
    plt.xlim((0.1, 100.0))
    plt.ylim((15.0, 70.0))
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.grid()
    plt.show()

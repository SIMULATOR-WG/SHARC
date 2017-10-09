# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 10:31:00 2017

@author: Andre Noll Barreto
"""

from sharc.propagation.propagation import Propagation
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


class PropagationBuildingEntryLoss(Propagation):
    """
    Implements the building entry loss according to ITU-R P.2109-0 (Prediction of Building Entry Loss)
    """

    def __init__(self):
        super().__init__()

    def get_loss(self, frequency_MHz, elevation, prob="random",
                 building_class="traditional", test = False) -> np.array:
        """
        Calculates building loss

        Parameters
        ----------
            frequency_MHz (np.array) : carrier frequencies (MHz)
            elevation (np.array) : apparent elevation angles
            prob (np.array / string) : the probability with which the loss is not exceeded;
                                    if "random", then different values are chosen for each user
            building_class (string) : type of construction material, "traditional" or "thermally-efficient"
            test (bool): True if only mu_1 is returned, for testing purposes (default False)

        Returns
        -------
            array with building loss values with dimensions of elevation
        """

        f_GHz = frequency_MHz / 1000

        if isinstance(prob, str) and prob.lower() == "random":
            prob = np.random.random(elevation.shape)

        if building_class == "traditional":
            r = 12.64
            s = 3.72
            t = .96
            u = 9.6
            v = 2.
            w = 9.1
            x = -3.
            y = 4.5
            z = -2.
        elif building_class == "thermally-efficient":
            r = 28.19
            s = -3.
            t = 8.48
            u = 13.5
            v = 3.8
            w = 27.8
            x = -2.9
            y = 9.4
            z = -2.1
        else:
            error_message = "building_class not supported"
            raise ValueError(error_message)
        c_dB = -3.
        hor_loss = r + s * np.log10(f_GHz) + t * (np.log10(f_GHz) ** 2)
        angle_correction = .212 * np.abs(elevation)
        mu_1 = hor_loss + angle_correction
        mu_2 = w + x * np.log10(f_GHz)
        sigma_1 = u + v * np.log10(f_GHz)
        sigma_2 = y + z * np.log10(f_GHz)

        a_dB = norm.ppf(prob) * sigma_1 + mu_1
        b_dB = norm.ppf(prob) * sigma_2 + mu_2

        a_lin = 10 ** (a_dB / 10)
        b_lin = 10 ** (b_dB / 10)
        c_lin = 10 ** (c_dB / 10)

        if test:
            b_lin = 0
            c_lin = 0

        loss = 10 * np.log10(a_lin + b_lin + c_lin)

        return loss

if __name__ == '__main__':

    entry_loss = PropagationBuildingEntryLoss()

    freq_GHz_log = np.arange(-1,2.1,.1)
    freq_GHz = 10 ** freq_GHz_log
    freq_MHz = freq_GHz * 1000

    # Plot median BLE mu_1, for comparison with ITU-R P2109-0
    plt.figure()
    median_loss_traditional = entry_loss.get_loss( freq_MHz, 0, prob=.5,
                                                   building_class="traditional", test=True)
    median_loss_therm_eff = entry_loss.get_loss(freq_MHz, 0, prob=.5,
                                                building_class="thermally-efficient", test=True)
    plt.semilogx(freq_GHz, median_loss_traditional, '-', label="traditional, 0deg")
    plt.semilogx(freq_GHz, median_loss_therm_eff, '--', label="thermally-efficient, 0deg")

    plt.legend(title="Building Type, elevation")
    plt.grid()
    plt.xlabel("frequency(GHz)")
    plt.ylabel("median loss (dB)")
    plt.title("Median Building Entry Loss (mu_1) - horizontal entry")



    # Plot median loss at different angles,
    # 0 degrees
    plt.figure()
    median_loss_traditional = entry_loss.get_loss( freq_MHz, 0, prob=.5,
                                                   building_class="traditional")
    median_loss_therm_eff = entry_loss.get_loss(freq_MHz, 0, prob=.5,
                                                building_class="thermally-efficient")

    plt.semilogx(freq_GHz, median_loss_traditional, '-', label="traditional, 0deg")
    plt.semilogx(freq_GHz, median_loss_therm_eff, '--', label="thermally-efficient, 0deg")

    # 45 degrees
    median_loss_traditional = entry_loss.get_loss( freq_MHz, 45, prob=.5,
                                                   building_class="traditional")
    median_loss_therm_eff = entry_loss.get_loss(freq_MHz, 45, prob=.5,
                                                building_class="thermally-efficient")

    plt.semilogx(freq_GHz, median_loss_traditional, '-', label="traditional, 45deg")
    plt.semilogx(freq_GHz, median_loss_therm_eff, '--', label="thermally-efficient, 45deg")

    # 90 deg
    median_loss_traditional = entry_loss.get_loss( freq_MHz, 90, prob=.5,
                                                   building_class="traditional")
    median_loss_therm_eff = entry_loss.get_loss(freq_MHz, 90, prob=.5,
                                                building_class="thermally-efficient")

    plt.semilogx(freq_GHz, median_loss_traditional, '-', label="traditional, 90deg")
    plt.semilogx(freq_GHz, median_loss_therm_eff, '--', label="thermally-efficient, 90deg")

    plt.legend(title="Building Type, elevation")
    plt.grid()
    plt.xlabel("frequency(GHz)")
    plt.ylabel("median loss (dB)")
    plt.title("Median Building Entry Loss - horizontal entry")


    plt.show()




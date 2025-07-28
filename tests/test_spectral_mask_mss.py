# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:56:10 2017

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.mask.spectral_mask_mss import SpectralMaskMSS


class SpectalMaskMSSTest(unittest.TestCase):
    def test_power_calc(self):
        # Test 1
        p_tx_density = -30  # dBW / Hz
        freq = 2190
        band = 1
        p_tx = p_tx_density + 10 * np.log10(band * 1e6)

        # reference bw is 4khz below 15GHz center freq
        p_tx_over_4khz = p_tx_density + 10 * np.log10(4e3)

        # dBm/MHz
        spurious_emissions = -30

        mask = SpectralMaskMSS(
            freq, band, spurious_emissions,
        )
        mask.set_mask(p_tx)

        N = len(mask.delta_f_lim)
        # N = 2

        should_eq = np.zeros(2 * N)
        eq = np.zeros(2 * N)
        spurious_start = 2 * band
        for i in range(N):
            f_offset = band / 2 + i * 4e-3

            F = (f_offset - band / 2) / band * 100

            should_eq[2 * i] = p_tx_over_4khz - 40 * np.log10(F / 50 + 1) + 30
            eq[2 * i] = mask.power_calc(freq + f_offset + 0.5 * 4e-3, 4e-3)

            should_eq[2 * i + 1] = should_eq[2 * i]
            eq[2 * i + 1] = mask.power_calc(freq - f_offset - 0.5 * 4e-3, 4e-3)

        # substitute last should eq with spurious emissions instead of formula
        should_eq[-1] = spurious_emissions + 10 * np.log10(4e-3)
        should_eq[-2] = spurious_emissions + 10 * np.log10(4e-3)

        npt.assert_almost_equal(should_eq, eq)

        npt.assert_equal(
            -np.inf,
            mask.power_calc(
                freq, band,
            ),
        )

        fll = mask.power_calc(freq + 1.5 * band, 2 * band)

        # this test only passes when considering 0 decimal places...
        # there is a noticeable difference in result from continous integration (by hand)
        # and the discrete integration implemented on SHARC. Could be mitigated
        # by implementing the mask differently
        self.assertAlmostEqual(fll, 10 * np.log10(165.33 * 1000), places=0)


if __name__ == '__main__':
    unittest.main()

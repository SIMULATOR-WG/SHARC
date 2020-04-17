from sharc.antenna.antenna_tactical_radio_relays import AntennaTacticalRadioRelays

import unittest
from numpy import testing as npt
import numpy as np


class TestAntennaTacticalRadioRelays(unittest.TestCase):
    def setUp(self) -> None:
        self.antenna = AntennaTacticalRadioRelays(None)

    def test_calculate_gain(self):
        # Test 1
        psi = np.array([0, 10, 19, 22, 89, 91, 120])
        gain = self.antenna.calculate_gain(off_axis_angle_vec=psi)
        expected_gain = np.array([20, 20, 20, 11, 11, 2, 2])
        npt.assert_equal(gain, expected_gain)


if __name__ == '__main__':
    unittest.main()

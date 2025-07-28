import unittest
import numpy as np
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.support.enumerations import StationType


class TestPropagationClutterLoss(unittest.TestCase):
    def setUp(self):
        self.clutter_loss = PropagationClutterLoss(np.random.RandomState(42))

    def test_spatial_clutter_loss(self):
        frequency = np.array([27000])  # MHz
        elevation = np.array([0, 45, 90])
        loc_percentage = np.array([0.1, 0.5, 0.9])
        distance = np.array([1000])  # meters, dummy value

        loss = self.clutter_loss.get_loss(
            distance=distance,
            frequency=frequency,
            elevation=elevation,
            loc_percentage=loc_percentage,
            station_type=StationType.FSS_SS,
        )

        # Check the shape of the output
        self.assertEqual(loss.shape, (3,))

        # Check if loss decreases with increasing elevation
        self.assertTrue(loss[0] >= loss[1] >= loss[2])

    def test_terrestrial_clutter_loss(self):
        frequency = np.array([2000, 6000])  # MHz
        distance = np.array([500, 2000])  # meters
        # Using a single value for location percentage
        loc_percentage = np.array([0.5])

        loss = self.clutter_loss.get_loss(
            frequency=frequency,
            distance=distance,
            loc_percentage=loc_percentage,
            station_type=StationType.IMT_BS,
        )

        self.assertEqual(loss.shape, (2,))

        self.assertTrue(loss[1] >= loss[0])

    def test_random_loc_percentage(self):
        frequency = np.array([4000])  # MHz
        distance = np.array([1000])  # meters

        loss = self.clutter_loss.get_loss(
            frequency=frequency,
            distance=distance,
            loc_percentage="RANDOM",
            station_type=StationType.IMT_UE,
        )

        self.assertTrue(0 <= loss <= 100)


if __name__ == '__main__':
    unittest.main()

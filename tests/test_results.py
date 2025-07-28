import unittest

from sharc.results import Results


class StationTest(unittest.TestCase):
    def setUp(self):
        self.results = Results().prepare_to_write(
            None,
            True,
            output_dir="output",
            output_dir_prefix="out",
        )

    def test_flush_to_and_load_from_file(self):
        arr1 = [1., 2., 3., 4., 5., 6., 7., 8., 9., 100.]
        self.results.imt_coupling_loss.extend(arr1)
        self.results.imt_bs_antenna_gain.extend(arr1)
        self.assertGreater(len(self.results.imt_coupling_loss), 0)
        self.assertGreater(len(self.results.imt_bs_antenna_gain), 0)
        # Results should flush
        self.results.write_files(1)
        # check that no results are left in arr
        self.assertEqual(len(self.results.imt_coupling_loss), 0)
        self.assertEqual(len(self.results.imt_bs_antenna_gain), 0)

        arr2 = [101., 102., 103., 104., 105., 106., 107., 108., 109.]
        self.results.imt_coupling_loss.extend(arr2)
        self.results.imt_bs_antenna_gain.extend(arr2)
        self.assertGreater(len(self.results.imt_coupling_loss), 0)
        self.assertGreater(len(self.results.imt_bs_antenna_gain), 0)
        self.results.write_files(2)
        # check that no results are left in arr
        self.assertEqual(len(self.results.imt_coupling_loss), 0)
        self.assertEqual(len(self.results.imt_bs_antenna_gain), 0)

        results_recuperated_from_file = Results().load_from_dir(self.results.output_directory)

        results_arr = arr1
        results_arr.extend(arr2)

        self.assertEqual(results_recuperated_from_file.imt_coupling_loss, results_arr)
        self.assertEqual(results_recuperated_from_file.imt_bs_antenna_gain, results_arr)

        results_recuperated_from_file = Results().load_from_dir(
            self.results.output_directory, only_samples=["imt_bs_antenna_gain"],
        )

        self.assertEqual(results_recuperated_from_file.imt_coupling_loss, [])
        self.assertEqual(results_recuperated_from_file.imt_bs_antenna_gain, results_arr)

    def test_get_most_recent_dirs(self):
        dir_2024_01_01_04 = "caminho_abs/prefixo_2024-01-01_04"
        dir_2024_01_01_10 = "caminho_abs/prefixo_2024-01-01_10"

        dir_2024_01_02_01 = "caminho_abs/prefixo_2024-01-02_01"
        dir_2024_10_01_01 = "caminho_abs/prefixo_2024-10-01_01"
        another_dir = "caminho_abs/prefixo2_2024-10-01_01"

        dirs = self.results.get_most_recent_outputs_for_each_prefix([
            dir_2024_01_01_04,
            dir_2024_01_01_10,
            dir_2024_01_02_01,
            dir_2024_10_01_01,
            another_dir,
        ])

        self.assertEqual(len(dirs), 2)

        self.assertIn(dir_2024_10_01_01, dirs)
        self.assertIn(another_dir, dirs)

        dirs = self.results.get_most_recent_outputs_for_each_prefix([
            dir_2024_01_01_04,
            dir_2024_01_01_10,
            dir_2024_01_02_01,
        ])

        self.assertEqual(len(dirs), 1)

        self.assertIn(dir_2024_01_02_01, dirs)

        dirs = self.results.get_most_recent_outputs_for_each_prefix([
            dir_2024_01_01_04,
            dir_2024_01_01_10,
        ])

        self.assertEqual(len(dirs), 1)

        self.assertIn(dir_2024_01_01_10, dirs)

    def test_get_prefix_date_and_id(self):
        dir_2024_01_01_04 = "caminho_abs/prefixo_2024-01-01_04"
        prefix, date, id = Results.get_prefix_date_and_id(dir_2024_01_01_04)
        self.assertEqual(prefix, "caminho_abs/prefixo_")
        self.assertEqual(date, "2024-01-01")
        self.assertEqual(id, "04")


if __name__ == '__main__':
    unittest.main()

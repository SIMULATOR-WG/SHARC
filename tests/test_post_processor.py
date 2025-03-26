
import unittest

from sharc.results import Results
from sharc.post_processor import PostProcessor


class StationTest(unittest.TestCase):
    def setUp(self):
        self.post_processor = PostProcessor()
        # We just prepare to write beacause Results class is not fully initialized
        # before preparing to read or loading from previous results
        self.results = Results().prepare_to_write(
            None,
            True
        )

    def test_generate_and_add_cdf_plots_from_results(self):
        self.results.imt_coupling_loss.extend([0, 1, 2, 3, 4, 5])
        self.results.imt_dl_inr.extend([0, 1, 2, 3, 4, 5])

        trace_legend = "any legendd. Lorem ipsum"
        self.post_processor.add_plot_legend_pattern(dir_name_contains="output", legend=trace_legend)

        self.post_processor.add_plots(
            self.post_processor.generate_cdf_plots_from_results(
                [self.results]
            )
        )

        self.assertEqual(len(self.post_processor.plots), 2)
        self.assertEqual(self.post_processor.plots[0].data[0].name, trace_legend)


if __name__ == '__main__':
    unittest.main()

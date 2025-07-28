import os
from pathlib import Path
from sharc.results import Results
from sharc.post_processor import PostProcessor

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
for i in range(0, 70, 10):
    post_processor.add_plot_legend_pattern(
        dir_name_contains="_lat_" + str(i) + "_deg",
        legend="latitude=" + str(i) + "deg")

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"), only_latest=True)
# ^: typing.List[Results]

post_processor.add_results(many_results)

plots = post_processor.generate_ccdf_plots_from_results(
    many_results
)

post_processor.add_plots(plots)

# Add a protection criteria line:
protection_criteria = -6
post_processor\
    .get_plot_by_results_attribute_name("imt_dl_inr", plot_type='ccdf')\
    .add_vline(protection_criteria, line_dash="dash")

# Show a single plot:
post_processor\
    .get_plot_by_results_attribute_name("imt_system_antenna_gain", plot_type='ccdf')\
    .show()

post_processor\
    .get_plot_by_results_attribute_name("system_imt_antenna_gain", plot_type='ccdf')\
    .show()

post_processor\
    .get_plot_by_results_attribute_name("sys_to_imt_coupling_loss", plot_type='ccdf')\
    .show()

post_processor\
    .get_plot_by_results_attribute_name("imt_system_path_loss", plot_type='ccdf')\
    .show()

post_processor\
    .get_plot_by_results_attribute_name("imt_dl_inr", plot_type='ccdf')\
    .show()

for result in many_results:
    # This generates the mean, median, variance, etc
    stats = PostProcessor.generate_statistics(
        result=result
    ).write_to_results_dir()

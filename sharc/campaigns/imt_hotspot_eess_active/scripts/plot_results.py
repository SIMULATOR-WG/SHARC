import os
from pathlib import Path
from sharc.results import Results
import plotly.graph_objects as go
from sharc.post_processor import PostProcessor
import pandas

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
post_processor\
    .add_plot_legend_pattern(
        dir_name_contains="beam_small_dl",
        legend="Small Beam DL"
    ).add_plot_legend_pattern(
        dir_name_contains="beam_small_ul",
        legend="Small Beam UL"
    )

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"), only_latest=True)

post_processor.add_results(many_results)

plots = post_processor.generate_cdf_plots_from_results(
    many_results
)

post_processor.add_plots(plots)

# This function aggregates IMT downlink and uplink
aggregated_results = PostProcessor.aggregate_results(
    dl_samples=post_processor.get_results_by_output_dir("_dl").system_inr,
    ul_samples=post_processor.get_results_by_output_dir("_ul").system_inr,
    ul_tdd_factor=0.25,
    # SF is not exactly 1, but approx
    n_bs_sim=1,
    n_bs_actual=1
)

# Add a protection criteria line:
# protection_criteria = int

# post_processor\
#     .get_plot_by_results_attribute_name("system_dl_interf_power")\
#     .add_vline(protection_criteria, line_dash="dash")

# Show a single plot:
relevant = post_processor\
    .get_plot_by_results_attribute_name("system_inr")

aggr_x, aggr_y = PostProcessor.cdf_from(aggregated_results)

relevant.add_trace(
    go.Scatter(x=aggr_x, y=aggr_y, mode='lines', name='Aggregate interference',),
)

compare_to = pandas.read_csv(
    os.path.join(campaign_base_dir, "comparison", "contribution_20.csv"),
    skiprows=1
)

comp_x, comp_y = (compare_to.iloc[:, 0], compare_to.iloc[:, 1])

relevant.add_trace(
    go.Scatter(x=comp_x, y=comp_y, mode='lines', name='Contrib 20 INR',),
)

relevant.show()

for result in many_results:
    # This generates the mean, median, variance, etc
    stats = PostProcessor.generate_statistics(
        result=result
    ).write_to_results_dir()

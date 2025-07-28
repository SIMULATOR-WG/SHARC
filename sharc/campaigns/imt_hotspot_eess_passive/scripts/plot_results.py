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
        dir_name_contains="1_cluster_DL_alternative",
        legend="Alternative Downlink"
    ).add_plot_legend_pattern(
        dir_name_contains="1_cluster_DL",
        legend="Downlink"
    ).add_plot_legend_pattern(
        dir_name_contains="1_cluster_UL",
        legend="Uplink"
    )

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

many_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output"), only_latest=True)

post_processor.add_results(many_results)

plots = post_processor.generate_cdf_plots_from_results(
    many_results
)

post_processor.add_plots(plots)

uplink_interf_samples = post_processor.get_results_by_output_dir("_UL").system_ul_interf_power

# This function aggregates IMT downlink and uplink
aggregated_results = PostProcessor.aggregate_results(
    dl_samples=post_processor.get_results_by_output_dir("_DL_alternative").system_dl_interf_power,
    ul_samples=uplink_interf_samples,
    ul_tdd_factor=0.25,
    # SF is not exactly 3, but approx
    n_bs_sim=1,
    n_bs_actual=3
)

# Add a protection criteria line:
# protection_criteria = int

# post_processor\
#     .get_plot_by_results_attribute_name("system_dl_interf_power")\
#     .add_vline(protection_criteria, line_dash="dash")

# Show a single plot:
relevant = post_processor\
    .get_plot_by_results_attribute_name("system_dl_interf_power")

# Title of CDF updated because ul interf power will be included
relevant.update_layout(
    title='CDF Plot for Interference',
)

aggr_x, aggr_y = PostProcessor.cdf_from(uplink_interf_samples)

relevant.add_trace(
    go.Scatter(x=aggr_x, y=aggr_y, mode='lines', name='Uplink',),
)

aggr_x, aggr_y = PostProcessor.cdf_from(aggregated_results)

relevant.add_trace(
    go.Scatter(x=aggr_x, y=aggr_y, mode='lines', name='Aggregate interference',),
)

# TODO: put some more stuff into PostProcessor if ends up being really used
compare_to = pandas.read_csv(
    os.path.join(campaign_base_dir, "comparison", "Fig. 8 EESS (Passive) Sensor.csv"),
    skiprows=1
)

comp_x, comp_y = (compare_to.iloc[:, 0], compare_to.iloc[:, 1])
# inverting given chart from P of I > x to P of I < x
comp_y = 1 - comp_y
# converting dB to dBm
comp_x = comp_x + 30

relevant.add_trace(
    go.Scatter(x=comp_x, y=comp_y, mode='lines', name='Fig. 8 EESS (Passive) Sensor',),
)

compare_to = pandas.read_csv(
    os.path.join(campaign_base_dir, "comparison", "Fig. 15 (IMT Uplink) EESS (Passive) Sensor.csv"),
    skiprows=1
)

comp_x, comp_y = (compare_to.iloc[:, 0], compare_to.iloc[:, 1])
# inverting given chart from P of I > x to P of I < x
comp_y = 1 - comp_y
# converting dB to dBm
comp_x = comp_x + 30

relevant.add_trace(
    go.Scatter(x=comp_x, y=comp_y, mode='lines', name='Fig. 15 (IMT Uplink) EESS (Passive) Sensor',),
)

relevant.show()

post_processor\
    .get_plot_by_results_attribute_name("system_ul_interf_power").show()

# for result in many_results:
#     # This generates the mean, median, variance, etc
#     stats = PostProcessor.generate_statistics(
#         result=result
#     ).write_to_results_dir()

aggregated_res_statistics = PostProcessor.generate_sample_statistics(
    "Aggregate Results Statistics",
    aggregated_results
)

print("\n###########################")

print(aggregated_res_statistics)

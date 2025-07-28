"""Script to process and plot results for MSS D2D to IMT cross-border campaign."""

import os
from pathlib import Path
from sharc.results import Results
# import plotly.graph_objects as go
from sharc.post_processor import PostProcessor

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file

prefixes = ["0km", "157.9km", "213.4km", "268.9km", "324.4km", "379.9km", "border"]
for prefix in prefixes:
    if prefix == "border":
        km = "0km"
    else:
        km = prefix
    post_processor\
        .add_plot_legend_pattern(
            dir_name_contains=f"{prefix}_base",
            legend=f"19 sectors ({km})"
        ).add_plot_legend_pattern(
            dir_name_contains=f"{prefix}_activate_random_beam_5p",
            legend=f"19 sectors, load=1/19 ({km})"
        ).add_plot_legend_pattern(
            dir_name_contains=f"{prefix}_activate_random_beam_30p",
            legend=f"19 sectors, load=30% ({km})"
        ).add_plot_legend_pattern(
            dir_name_contains=f"{prefix}_random_pointing_1beam",
            legend=f"1 sector random pointn ({km})"
        )

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())

results_dl = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output_base_dl"), only_latest=True)
results_ul = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output_base_ul"), only_latest=True)
# ^: typing.List[Results]
all_results = [*results_ul, *results_dl]

post_processor.add_results(all_results)

styles = ["solid", "dot", "dash", "longdash", "dashdot", "longdashdot"]


def linestyle_getter(result: Results):
    """
    Returns a line style string based on the prefix found in the result's output directory.
    """
    for i in range(len(prefixes)):
        if prefixes[i] in result.output_directory:
            return styles[i]
    return "solid"


post_processor.add_results_linestyle_getter(linestyle_getter)

plots = post_processor.generate_cdf_plots_from_results(
    all_results
)

post_processor.add_plots(plots)

# Add a protection criteria line:
protection_criteria = -6
post_processor\
    .get_plot_by_results_attribute_name("imt_dl_inr")\
    .add_vline(protection_criteria, line_dash="dash")

post_processor\
    .get_plot_by_results_attribute_name("imt_ul_inr")\
    .add_vline(protection_criteria, line_dash="dash")

# Add a protection criteria line:
pfd_protection_criteria = -109
post_processor\
    .get_plot_by_results_attribute_name("imt_dl_pfd_external_aggregated")\
    .add_vline(pfd_protection_criteria, line_dash="dash")


attributes_to_plot = [
    # "imt_system_antenna_gain",
    # "system_imt_antenna_gain",
    # "sys_to_imt_coupling_loss",
    # "imt_system_path_loss",
    "imt_dl_pfd_external",
    "imt_dl_pfd_external_aggregated",
    "imt_dl_inr",
    "imt_ul_inr",
]

for attr in attributes_to_plot:
    post_processor.get_plot_by_results_attribute_name(attr).show()

# Ensure the "htmls" directory exists relative to the script directory
# htmls_dir = Path(__file__).parent / "htmls"
# htmls_dir.mkdir(exist_ok=True)
# for attr in attributes_to_plot:
#     post_processor\
#         .get_plot_by_results_attribute_name(attr)\
#         .write_html(htmls_dir / f"{attr}.html")

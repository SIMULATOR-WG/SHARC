from sharc.run_multiple_campaigns_mut_thread import run_campaign_re

# Set the campaign name
# The name of the campaign to run. This should match the name of the campaign directory.
name_campaign = "imt_ntn_to_imt_tn_co_channel"

# Run the campaigns
# This function will execute the campaign with the given name.
# It will look for the campaign directory under the specified name and start the necessary processes.
regex_pattern = r'^parameters_imt_ntn_to_imt_tn_co_channel_sep_.*_km.yaml'
run_campaign_re("imt_ntn_to_imt_tn_co_channel", regex_pattern)
print("Executing campaign with regex pattern:", regex_pattern)

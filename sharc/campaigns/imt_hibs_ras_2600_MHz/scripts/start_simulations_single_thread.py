from sharc.run_multiple_campaigns import run_campaign

# Set the campaign name
# The name of the campaign to run. This should match the name of the campaign directory.
name_campaign = "imt_hibs_ras_2600_MHz"

# Run the campaign in single-thread mode
# This function will execute the campaign with the given name in a single-threaded manner.
# It will look for the campaign directory under the specified name and start the necessary processes.
run_campaign(name_campaign)

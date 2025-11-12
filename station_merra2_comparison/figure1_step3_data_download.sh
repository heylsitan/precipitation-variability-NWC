#!/bin/bash

# Configure input and output paths
input_dir=""
output_dir=""  # Please modify to actual output path

# Automatically create output directory
mkdir -p ${output_dir}

# Define start and end years
start_year=1982
end_year=2020

# Loop through each year
for year in $(seq $start_year $end_year); do
    # Define input file matching pattern
    input_files=""
    
    # Define output filename
    output_file=""
    
    # Check if matching files exist
    if ls ${input_files} 1> /dev/null 2>&1; then
        echo "Merging data for year ${year}..."
        # Use CDO to merge files (with memory optimization parameters)
        cdo  mergetime ${input_files} ${output_file}
        echo "Year ${year} merge completed -> ${output_file}"
    else
        echo "[Warning] No data files for year ${year}, skipping"
    fi
done
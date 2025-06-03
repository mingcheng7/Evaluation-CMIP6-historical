#!/bin/bash

# Output directory
work_dir="wk"
mkdir -p $work_dir
output_dir="processed_files"
mkdir -p $output_dir

# Run loops
for var in $(ls *.nc | awk -F_ '{print $1}' | sort | uniq); do
        for model in $(ls ${var}_*.nc | awk -F_ '{print $3}' | sort | uniq); do
                echo "Processing variable: $var, model: $model"

                files=$(ls ${var}_*_${model}_*.nc)

                merged_file="${work_dir}/${var}_${model}_merged.nc"
                cdo mergetime $files $merged_file

                time_clipped_file="${work_dir}/${var}_${model}_time_clipped.nc"
                cdo seldate,2000-01-01,2014-12-31 $merged_file $time_clipped_file

                regridded_file="${work_dir}/${var}_${model}_regridded.nc"
                cdo remapbil,r360x180 $time_clipped_file $regridded_file

                final_file="${work_dir}/${var}_${model}_final.nc"
                cdo sellonlatbox,-360,360,-90,-30 $regridded_file $final_file

                new_file="${out_dir}/${var}_Omon_${model}_historical_r1i1p1f1_rg_so_200001-201412.nc"
                cp $final_file $new_file

                echo "Finished processing $newfile"
        done
done
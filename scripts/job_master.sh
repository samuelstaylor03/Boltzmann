#!/bin/bash

# Define variables
r_start=1
r_end=5
temp=300
cur_runs_dir="y_0.5_runs"
copy_dir="C2H2_copy"
control_file_name="control_boltzmann.inp"
job_file_name="job_script.sh"
boltzmann_script="./boltzmann"
job_name_pre_r_val="C2H2_y_0.5_r"

# Function to copy directories, update files, and submit jobs
process_and_submit() {
    for (( r_val=r_start; r_val<=r_end; r_val++ )); do
        # New directory name
        new_dir="r${r_val}"

        # Check if the directory already exists
        if [ -d "${new_dir}" ]; then
            echo "Directory ${new_dir} already exists, skipping copy"
            continue
        fi

        # Copy the directory
        if ! cp -r "${copy_dir}" "${new_dir}"; then
            echo "Error copying ${copy_dir} to ${new_dir}"
            exit 1
        fi
        echo "Copied ${copy_dir} to ${new_dir}"

        # Update the velocity_output_path in the control_boltzmann.inp file
        sed -i "s|^velocity_output_path=.*|velocity_output_path=/scratch/group/p.phy240167.000/collision/C2H2/${cur_runs_dir}/r${r_val}/|" "${control_file_name}"
        echo "Updated velocity_output_path in ${control_file_name}"

        # Run the boltzmann script with the current r_val
        if ! "${boltzmann_script}" "${temp}" "${r_val}"; then
            echo "Error running boltzmann script for r${r_val}"
            exit 1
        fi
        echo "Ran boltzmann script for r${r_val}"

        # Change to the new directory
        cd "${new_dir}" || { echo "Failed to change to directory ${new_dir}"; exit 1; }

        # Update the job name in the job_script.sh file
        sed -i "s/^#SBATCH --job-name=.*/#SBATCH --job-name=${job_name_pre_r_val}${r_val}/" "${job_file_name}"
        echo "Updated job name in ${job_file_name}"

        # Submit the job script
        if ! sbatch "${job_file_name}"; then
            echo "Error submitting job for r${r_val}"
            exit 1
        fi
        echo "Submitted job for r${r_val}"

        # Change back to the original directory
        cd .. || { echo "Failed to return to the parent directory"; exit 1; }
    done

    echo "All jobs processed and submitted."
}

# Main execution starts here
process_and_submit

exit 0

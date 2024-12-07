import os
import pkgutil
import io
from pathlib import Path

template_path = (Path(__file__) / '../templates' ).resolve()

# Tasks

# 1. Extract number of expected cells from metrics CSV per sample
# 2. Modify slurm file


def create_cellbender_slurm(
        conda_env,
        raw_h5,
        output_h5,
        expected_cells,
        slurm_job,
        slurm_time,
        slurm_mem,
        output_slurm_file        
):
    slurm_template = template_path / "cellbender_template.slurm"
    to_replace = [
        "USER_CONDA_ENV",
        "USER_RAW_MATRIX_FILE",
        "USER_OUTPUT_MATRIX_FILE",
        # "USER_EXPECTED_CELL_COUNT",
        "USER_JOB_NAME",
        "USER_TIME",
        "USER_MEM"
    ]

    value = [
        conda_env,
        str(raw_h5),
        str(output_h5),
        # str(expected_cells),
        slurm_job,
        slurm_time,
        str(slurm_mem) + "G"
    ]

    infile = open(slurm_template, 'r')
    outfile = open(output_slurm_file, 'w')

    for line in infile:
        for check, rep in zip(to_replace, value):
            line = line.replace(check, rep)
        outfile.write(line)
    infile.close()
    outfile.close()
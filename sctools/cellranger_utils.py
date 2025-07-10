import typer
import os
import pkgutil
import pandas as pd
import numpy as np
import subprocess as sp
import io
from pathlib import Path
from typing import Optional
from typing_extensions import Annotated

template_path = (Path(__file__) / '../templates' ).resolve()

def parse_metrics(
    sample,
    multi_output_dir,
):
    base_dir = multi_output_dir / "outs" / "per_sample_outs" / sample / "metrics_summary.csv"
    try:
        metrics_df = pd.read_csv(base_dir, skip_blank_lines=True, comment = "#")
    except FileNotFoundError as e:
        print(f"Error: {e}")

    cells = metrics_df.loc[metrics_df['Metric Name'] == 'Cells',"Metric Value"].iloc[0]
    cells = int(cells.replace(",",""))

    reads = metrics_df.loc[(metrics_df['Metric Name'] == 'Number of reads') & 
                      (metrics_df['Library Type'] == 'Gene Expression'), 'Metric Value'].iloc[0]
    reads = int(reads.replace(",",""))
    return(reads, cells)

def create_count_config(
        library,
        genome_reference,
        cells,
        feature_reference=None,
        vdj_reference=None,
        fastq_dir=None,
        BCR=None,
        TCR=None,
        antibody=None,
        outdir=None
):
    # Load the config template
    config = pd.read_csv(template_path / "sample_config_template.csv", header=None)

    # Drop rows for BCR, TCR, or antibody if values are None
    drop_rows = []

    if BCR is None or not Path(BCR).exists() or not any(Path(BCR).iterdir()):
        drop_rows += config[config.apply(lambda row: "USER_BCR" in ",".join(row.astype(str)), axis=1)].index.tolist()

    if TCR is None or not Path(TCR).exists() or not any(Path(TCR).iterdir()):
        drop_rows += config[config.apply(lambda row: "USER_TCR" in ",".join(row.astype(str)), axis=1)].index.tolist()

    if antibody is None or not Path(antibody).exists() or not any(Path(antibody).iterdir()):
        drop_rows += config[config.apply(lambda row: "USER_ANTIBODY" in ",".join(row.astype(str)), axis=1)].index.tolist()

    config.drop(index=drop_rows, inplace=True)

    # Replace placeholder values
    config.replace(
        to_replace=[
            "USER_REF", "USER_CELLS", "USER_FEATURE_REF", "USER_VDJ_REF",
            "USER_FASTQ_DIR", "USER_BCR_DIR", "USER_TCR_DIR", "USER_ANTIBODY_DIR",
            "USER_BCR_SAMPLE", "USER_TCR_SAMPLE", "USER_ANTIBODY_SAMPLE"
        ],
        value=[
            str(genome_reference),
            cells,
            str(feature_reference) if feature_reference else "",
            str(vdj_reference) if vdj_reference else "",
            str(fastq_dir),
            str(BCR) if BCR else "",
            str(TCR) if TCR else "",
            str(antibody) if antibody else "",
            library + "B",
            library + "T",
            library + "F"
        ],
        inplace=True
    )

    config.to_csv(outdir, header=False, index=False)


def count_cmd(
    cellranger_path: Annotated[str, typer.Option(help="Path to cellranger")] = "$GROUPDIR/$USER/tools/cellranger-9.0.0",
    config: Annotated[str, typer.Argument(help="Path to the config file")] = "configs/config.csv",
    sample: Annotated[str, typer.Argument(help="The sample to process")] = None,
    slurm_mode: Annotated[bool, typer.Option(help="Run cellranger using built-in SLURM mode")] = False,
    threads: Annotated[int, typer.Option(help="Number of CPU to use")] = 32,
    memory: Annotated[int, typer.Option(help="Memory to use in GB per CPU")] = 4,
):
    cellranger_cmd = [
        cellranger_path / "cellranger",
        "multi",
        "--id", sample,
        "--csv", config,
    ]
    if slurm_mode:
        cellranger_cmd = cellranger_cmd + ["--jobmode=slurm"]
    else:
        cellranger_cmd = cellranger_cmd + ["--localcores", str(threads), "--localmem", str(memory)]
        slurm_cmd = ["srun","--job-name="+sample+"_count",
                     "--ntasks=1", "--cpus-per-task="+str(threads),
                     "--mem="+str(memory)+"G", "--time=1:00:00",
                     "--output="+sample+"_count_%j.out",
                     "--error="+sample+"_count_%j.out"]
        cellranger_cmd = slurm_cmd+cellranger_cmd+["&"]
    
    print(" ".join([str(x) for x in cellranger_cmd]))
    sp.Popen(cellranger_cmd)


def get_dir_size(path='.'):
    total = 0
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_file():
                total += entry.stat().st_size
            elif entry.is_dir():
                total += get_dir_size(entry.path)
    return total
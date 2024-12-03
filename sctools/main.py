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

app = typer.Typer()
template_path = (Path(__file__) / '../templates' ).resolve()

# Python CLI

# You must first open an interactive node (minimal resources) on slurm
# and then run the CLI which will automatically execute each step of the pipeline
# using SLURM on the HPC

@app.command()
def demux(
    cellranger_path: Annotated[Optional[Path], typer.Option(help="Path to cellranger")] = "$GROUPDIR/$USER/tools/cellranger-9.0.0",
    config: Annotated[Optional[Path], typer.Argument(help="Path to the demultiplexing config file")] = "configs/config.csv",
    library: Annotated[str, typer.Argument(help="The library to process")] = None,
    slurm_mode: Annotated[bool, typer.Option(help="Run cellranger using built-in SLURM mode")] = False,
    threads: Annotated[int, typer.Option(help="Number of CPU to use")] = 32,
    memory: Annotated[int, typer.Option(help="Memory to use in GB per CPU")] = 4,
    demux_output_suffix: Annotated[str, typer.Option(help="Suffix to add to the demux output directory")] = "_DEMUX",
    bamtofastq_dir: Annotated[str, typer.Option(help="Directory to save bamtofastq output files")] = "bamtofastq",
    genome_reference: Annotated[str, typer.Option(help="Path to the genome reference")] = "$GROUPDIR/$USER/projects/scrna/human/references/refdata-gex-GRCh38-2024-A",
    feature_reference: Annotated[str, typer.Option(help="Path to the feature reference (e.g. Biolegend_feature.csv")] = "$GROUPDIR/$USER/projects/scrna/human/references/Biolegend_feature.csv",
    vdj_reference: Annotated[str, typer.Option(help="Path to the VDJ reference")] = "$GROUPDIR/$USER/projects/scrna/human/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
):
    """
    
    """
    # Ensure environment includes cellranger exec files
    cellranger_path = str(cellranger_path)
    
    # Read the config file
    config_df = pd.read_csv(config, skip_blank_lines=True, comment = "#")

    # Create cellranger multi command
    crmulti_cmd = [cellranger_path+"/cellranger","multi",
                   "--id", library+"_DEMUX",
                   "--csv", config]
    if slurm_mode:
        crmulti_cmd.append("--jobmode=slurm")
    else:
        slurm_cmd = ["srun","--job-name="+library+"_demux",
                     "--ntasks=1", "--cpus-per-task="+str(threads),
                     "--mem="+str(memory)+"G", "--time=1:00:00",
                     "--output="+library+"_demux_%j.out",
                     "--error="+library+"_demux_%j.out"]
        crmulti_cmd = slurm_cmd+crmulti_cmd
    
    # Run cellranger multi
    sp.run(crmulti_cmd)

    # Record metrics for each sample and run bamtofastq
    Path("bamtofastq").mkdir(parents=True, exist_ok=True)

    sample_index = list(config_df['[gene-expression]']).index("sample_id")
    samples = list(config_df.iloc[sample_index+1:,0])
    sample_cells = dict()
    sample_reads = dict()
    sample_bamtofastq_cmd = []
    for sample in samples:
        sample_metrics = parse_metrics(sample, multi_output_dir=library+"_DEMUX")
        sample_reads[sample] = sample_metrics[0]
        sample_cells[sample] = sample_metrics[1]

        sp.run(["echo","Creating job for sample: " + sample + "..."])

        # Round up reads to ensure you're counting everything
        reads_str = str(sample_reads[sample])
        first_two_digits = int(reads_str[:2])
        first_two_digits += 1
        rounded_reads = str(first_two_digits) + '0' * (len(reads_str) - 2)

        bamtofq_cmd = [
            "srun",
            "--cpus-per-task=8",
            "--mem=16G",
            "--time=0:30:00",
            f"--output=bamtofastq/{sample}.out",
            f"--error=bamtofastq/{sample}.err",
            f"--job-name={sample}_bamtofastq",
            f"{cellranger_path}/lib/bin/bamtofastq",  # The actual bamtofastq command starts here
            "--nthreads=8",  # Matches --cpus-per-task
            "--reads-per-fastq", str(rounded_reads),
            f"{library}_DEMUX/outs/per_sample_outs/{sample}/count/sample_alignments.bam",
            f"bamtofastq/{sample}"
        ]

        print("Running command:", " ".join(bamtofq_cmd))
        sp.Popen(bamtofq_cmd)
        # Create cellranger count config files
        create_count_config(
            library=library,
            sample_id=sample,
            cells=sample_cells[sample],
            genome_reference=genome_reference,
            feature_reference=feature_reference,
            vdj_reference=vdj_reference,
            GEX_prefix="bamtofastq",
            BCR="BCR",
            TCR="TCR",
            antibody="MC_AB",
            outdir="configs"
        )

    ### Final cellranger count
    for sample in samples:
        count(
            cellranger_path=cellranger_path,
            config="configs/"+sample+"_config.csv",
            sample=sample,
            threads=threads,
            memory=memory,
            slurm_mode=slurm_mode
        )
        


def parse_metrics(
    sample,
    multi_output_dir,
):
    base_dir = multi_output_dir + "/outs/per_sample_outs/" + sample + "/metrics_summary.csv"
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
        sample_id,
        genome_reference,
        cells,
        feature_reference = None,
        vdj_reference = None,
        GEX_prefix="bamtofastq",
        BCR=None,
        TCR=None,
        antibody=None,
        outdir="configs"
):        
    config = pd.read_csv(template_path / "sample_config_template.csv", 
                         skip_blank_lines=True, comment = "#", header=None)


    config.replace(
        to_replace=["USER_REF",
                    "USER_CELLS",
                    "USER_FEATURE_REF",
                    "USER_VDJ_REF",
                    "USER_BCL_DIR",
                    "USER_TCL_DIR",
                    "USER_AB_DIR"],
        value=[genome_reference,
               cells,
               feature_reference,
               vdj_reference,
               BCR,
               TCR,
               antibody],
        inplace=True
    )
    config.to_csv(outdir+"/"+sample_id+"_config.csv", header = None, index=False)

@app.command()
def count(
    cellranger_path: Annotated[str, typer.Option(help="Path to cellranger")] = "$GROUPDIR/$USER/tools/cellranger-9.0.0",
    config: Annotated[str, typer.Argument(help="Path to the config file")] = "configs/config.csv",
    sample: Annotated[str, typer.Argument(help="The sample to process")] = None,
    slurm_mode: Annotated[bool, typer.Option(help="Run cellranger using built-in SLURM mode")] = False,
    threads: Annotated[int, typer.Option(help="Number of CPU to use")] = 32,
    memory: Annotated[int, typer.Option(help="Memory to use in GB per CPU")] = 4,
):
    cellranger_cmd = [
        cellranger_path+"/cellranger",
        "multi",
        "--id", sample,
        "--csv", config,
        "--localcores", str(threads),
        "--localmem", str(memory),
    ]
    if slurm_mode:
        cellranger_cmd.append("--jobmode=slurm")
    else:
        slurm_cmd = ["srun","--job-name="+sample+"_count",
                     "--ntasks=1", "--cpus-per-task="+str(threads),
                     "--mem="+str(memory)+"G", "--time=1:00:00",
                     "--output="+sample+"_count_%j.out",
                     "--error="+sample+"_count_%j.out"]
        cellranger_cmd = slurm_cmd+cellranger_cmd
    
    sp.run(cellranger_cmd)

    

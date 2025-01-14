import sctools.cellranger_utils as cr_utils
import sctools.cellbender_utils as cb_utils
import typer
import os
import pkgutil
import pandas as pd
import numpy as np
import subprocess as sp
import io
import datetime
from pathlib import Path
from typing import Optional
from typing_extensions import Annotated

app = typer.Typer()

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
    for sample in samples:
        sample_metrics = cr_utils.parse_metrics(sample, multi_output_dir=Path(library+demux_output_suffix))
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
        



@app.command()
def crcount(
    library: Annotated[str,typer.Argument(help="Library name")] = None,
    cellranger_path: Annotated[
        Path, 
        typer.Option(
            exists = True,
            file_okay = False,
            dir_okay = True,
            help="Path to cellranger"
        ),
    ] = "$GROUPDIR/$USER/tools/cellranger-9.0.0",
    demux_path: Annotated[
        Path, 
        typer.Argument(
            exists = True,
            file_okay = False,
            dir_okay = True,
            help="Path to the demultiplexed output folder"
        ),
    ] = None,
    demux_config: Annotated[
        Path, 
        typer.Argument(
            exists = True,
            file_okay = True,
            dir_okay = False,
            help="Path to the demultiplexing config CSV file"
        ),
    ] = None,
    bamtofastq_dir: Annotated[
        Path, 
        typer.Argument(
            exists = True,
            file_okay = False,
            dir_okay = True,
            resolve_path = True,
            help="Path to the bamtofastq output folder"
        ),
    ] = Path("bamtofastq"),
    genome_reference: Annotated[
        Path, 
        typer.Option(
            exists = True,
            file_okay = False,
            dir_okay = True,
            help="Path to the genome reference"
        ),
    ] = "$GROUPDIR/$USER/projects/scrna/human/references/refdata-gex-GRCh38-2024-A",
    feature_reference: Annotated[
        Path, 
        typer.Option(
            exists = True,
            file_okay = True,
            dir_okay = False,
            help="Path to the feature reference (e.g. Biolegend_feature.csv"
        ),
    ] = "$GROUPDIR/$USER/projects/scrna/human/references/Biolegend_feature.csv",
    vdj_reference: Annotated[
        Path, 
        typer.Option(
            exists = True,
            file_okay = False,
            dir_okay = True,
            help="Path to the vdj reference"
        ),
    ] = Path("$GROUPDIR/$USER/projects/scrna/human/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"),
    BCR: Annotated[
        Path, 
        typer.Option(
            exists = True,
            file_okay = False,
            dir_okay = True,
            resolve_path = True,
            help="Path to the BCR FASTQs"
        ),
    ] = Path("BCR"),
    TCR: Annotated[
        Path, 
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            resolve_path = True,
            help="Path to the TCR FASTQs"
        ),
    ] = Path("TCR"),
    antibody: Annotated[
        Path, 
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            resolve_path = True,
            help="Path to the antibody FASTQs"
        ),
    ] = Path("MC_AB"),
    config_dir: Annotated[
        Path, 
        typer.Option(
            exists=True,
            file_okay=False,
            dir_okay=True,
            help="Path to the output folder for the config files"
        ),
    ] = Path("configs"),
    slurm_mode: Annotated[bool, typer.Option(help="Run cellranger using built-in SLURM mode")] = False,
    threads: Annotated[int, typer.Option(help="Number of CPU to use")] = 32,
    memory: Annotated[int, typer.Option(help="Memory to use in GB per CPU")] = 4
):
    # 1. Extract info to create config files
    demux_config_csv = pd.read_csv(demux_config, skip_blank_lines=True, comment = "#")
    sample_index = list(demux_config_csv['[gene-expression]']).index("sample_id")
    samples = list(demux_config_csv.iloc[sample_index+1:,0])        
    for sample in samples:
        sample_cells = cr_utils.parse_metrics(sample, multi_output_dir=demux_path)[1]

        subdirs = [x for x in (bamtofastq_dir / sample).iterdir() if x.is_dir()]        
        subdirs = dict(zip(subdirs, [cr_utils.get_dir_size(x) for x in subdirs]))
        GEX_fastq_dir = max(subdirs.items(), key=lambda item: item[1])[0]
        config_path = str(config_dir) + "/" + sample + "_config.csv"
        
        # 2. Create config files

        cr_utils.create_count_config(
            library=library,
            cells=sample_cells,
            genome_reference=genome_reference,
            feature_reference=feature_reference,
            vdj_reference=vdj_reference,
            fastq_dir=GEX_fastq_dir,
            BCR=BCR,
            TCR=TCR,
            antibody=antibody,
            outdir=config_path
        )

        cr_utils.count_cmd(
            cellranger_path=cellranger_path,
            config=config_path,
            sample=sample,
            slurm_mode=slurm_mode,
            threads=threads,
            memory=memory
        )

@app.command()
def cellbender(
    conda_env: Annotated[Optional[str], typer.Option(help = "Name of cellbender conda environment for your user")] = "cellbender",
    sample: Annotated[str, typer.Argument(help = "Sample name. This should correspond to the sample folder in the specified 'library_dir'")] = None,
    library_dir: Annotated[Path, typer.Argument(
        help = "Path to the raw multiplexed library that contains this sample",
        dir_okay = True, file_okay = False, exists = True
    )] = None,
    output_dir: Annotated[Path, typer.Option(
        help = "Path to the output parent directory for the cellbender results",
        dir_okay = True, exists = True
    )] = Path("/groups/ag2239_gp/projects/scrna/human/cellbender"),
    slurm_time: Annotated[Optional[str], typer.Option(help = "Time for the slurm job")] = "03:00:00",
    slurm_mem: Annotated[Optional[int], typer.Option(help = "Memory for the slurm job")] = 4,
    use_expected_cells = False,
):
    # 1. Extract info to create config files
    sample_dir = library_dir / sample
    if use_expected_cells:
        sample_cells = cr_utils.parse_metrics(sample, sample_dir)[1]
    else:
        sample_cells = None
    
    formatted_date = datetime.date.today().strftime("%Y_%m_%d")
    output_dir_full = output_dir / sample / (sample+"_"+formatted_date)
    output_dir_full.mkdir(parents=True, exist_ok=True)
    
    raw_h5_path = sample_dir.resolve() / 'outs' / 'multi' / 'count' / 'raw_feature_bc_matrix.h5'

    os.chdir(output_dir_full)
    cb_utils.create_cellbender_slurm(
        conda_env=conda_env,
        raw_h5 = raw_h5_path.resolve(),
        output_h5 = output_dir_full / 'cellbender_output_filtered.h5',
        expected_cells = sample_cells,
        slurm_job = sample + "_cellbender",
        slurm_time = slurm_time,
        slurm_mem = slurm_mem,
        output_slurm_file = "cellbender.slurm"
    )

    os.chdir(output_dir_full)
    sp.Popen(["sbatch", "cellbender.slurm"])

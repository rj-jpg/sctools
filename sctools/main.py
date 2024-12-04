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
        sample_metrics = parse_metrics(sample, multi_output_dir=library+demux_output_suffix)
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
        # create_count_config(
        #     library=library,
        #     sample_id=sample,
        #     cells=sample_cells[sample],
        #     genome_reference=genome_reference,
        #     feature_reference=feature_reference,
        #     vdj_reference=vdj_reference,
        #     GEX_prefix="bamtofastq",
        #     BCR="BCR",
        #     TCR="TCR",
        #     antibody="MC_AB",
        #     outdir="configs"
        # )

    ### Final cellranger count
    # for sample in samples:
    #     count(
    #         cellranger_path=cellranger_path,
    #         config="configs/"+sample+"_config.csv",
    #         sample=sample,
    #         threads=threads,
    #         memory=memory,
    #         slurm_mode=slurm_mode
    #     )
        


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
        sample_id,
        genome_reference,
        cells,
        feature_reference = None,
        vdj_reference = None,
        fastq_dir = None,
        BCR=None,
        TCR=None,
        antibody=None,
        outdir=None
):        
    config = pd.read_csv(template_path / "sample_config_template.csv", header=None)


    config.replace(
        to_replace=["USER_REF",
                    "USER_CELLS",
                    "USER_FEATURE_REF",
                    "USER_VDJ_REF",
                    "USER_FASTQ_DIR",
                    "USER_BCL_DIR",
                    "USER_TCL_DIR",
                    "USER_ANTIBODY_DIR"],
        value=[str(genome_reference),
               cells,
               str(feature_reference),
               str(vdj_reference),
               str(fastq_dir),
               str(BCR),
               str(TCR),
               str(antibody)],
        inplace=True
    )
    config.to_csv(outdir, header = False, index=False)

@app.command()
def crcount(
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
        sample_cells = parse_metrics(sample, multi_output_dir=demux_path)[1]

        subdirs = [x for x in (bamtofastq_dir / sample).iterdir() if x.is_dir()]        
        subdirs = dict(zip(subdirs, [get_dir_size(x) for x in subdirs]))
        GEX_fastq_dir = max(subdirs.items(), key=lambda item: item[1])[0]
        config_path = str(config_dir) + "/" + sample + "_config.csv"
        
        # 2. Create config files

        create_count_config(
            sample_id=sample,
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

        count_cmd(
            cellranger_path=cellranger_path,
            config=config_path,
            sample=sample,
            slurm_mode=slurm_mode,
            threads=threads,
            memory=memory
        )

    
    
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
        cellranger_cmd.append("--jobmode=slurm")
    else:
        cellranger_cmd = cellranger_cmd + ["--localcores", str(threads), "--localmem", str(memory)]
        slurm_cmd = ["srun","--job-name="+sample+"_count",
                     "--ntasks=1", "--cpus-per-task="+str(threads),
                     "--mem="+str(memory)+"G", "--time=1:00:00",
                     "--output="+sample+"_count_%j.out",
                     "--error="+sample+"_count_%j.out"]
        cellranger_cmd = slurm_cmd+cellranger_cmd
    
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

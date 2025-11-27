# FASTQ Processing Pipeline

A Nextflow pipeline for processing paired-end FASTQ files through quality control, trimming, and alignment.

## Pipeline Overview

This pipeline performs the following steps:

1. **FastQC (Raw)** - Quality control on raw FASTQ files
2. **fastp** - Adapter trimming and quality filtering
3. **FastQC (Trimmed)** - Quality control on trimmed reads
4. **BWA-MEM** - Alignment to reference genome
5. **BAM Sorting** - Sort and index aligned reads
6. **Alignment Statistics** - Calculate mapping statistics
7. **MultiQC** - Aggregate all QC reports into a single HTML report

## Requirements

- [Nextflow](https://www.nextflow.io/) (>=21.10.0)
- [Docker](https://www.docker.com/) (or Singularity)
- Reference genome FASTA file
- Paired-end FASTQ files

### Docker Images Used

The pipeline automatically pulls the following Docker containers:
- `biocontainers/fastqc:v0.11.9_cv8`
- `biocontainers/fastp:v0.23.4_cv1`
- `biocontainers/bwa:v0.7.17_cv1`
- `biocontainers/samtools:v1.19.2_cv1`
- `ewels/multiqc:v1.19`

## Quick Start

### 1. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
```

### 2. Prepare Your Input Files

#### Option A: Create Samplesheet Automatically

Use the provided Python script to automatically generate a samplesheet from a directory of FASTQ files:

```bash
# Install Python dependencies (none required - uses stdlib only)
pip install -r requirements.txt

# Generate samplesheet from a directory
python bin/create_samplesheet.py /path/to/fastq/directory

# This creates: /path/to/fastq/directory/<directory_name>_SampleSheet.csv
```

The script will:
- Recursively search for all FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz)
- Automatically detect paired-end reads (supports both `_R1/_R2` and `_1/_2` patterns)
- Extract sample names by removing standard Illumina suffixes
- Handle unpaired reads by placing them in the fastq1 column
- Generate a samplesheet compatible with the pipeline

**Usage examples:**
```bash
# Basic usage - creates samplesheet in the same directory
python bin/create_samplesheet.py /data/fastq_files

# Specify custom output location
python bin/create_samplesheet.py /data/fastq_files --output my_samples.csv

# Use relative paths instead of absolute paths
python bin/create_samplesheet.py /data/fastq_files --relative-paths

# Verbose output for debugging
python bin/create_samplesheet.py /data/fastq_files --verbose
```

#### Option B: Create Samplesheet Manually

Create a samplesheet CSV file with three columns: `sampleID`, `fastq1`, `fastq2`

```csv
sampleID,fastq1,fastq2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

### 3. Configure Parameters

Edit `params.yaml` to specify your reference genome and adjust parameters:

```yaml
reference: "/path/to/reference/genome.fa"
samplesheet: "samplesheet.csv"
outdir: "results"
```

### 4. Run the Pipeline

```bash
nextflow run main.nf -params-file params.yaml
```

## Usage

### Basic Usage

```bash
nextflow run main.nf \
    --samplesheet samplesheet.csv \
    --reference /path/to/reference.fa \
    --outdir results
```

### Using a Parameters File

```bash
nextflow run main.nf -params-file params.yaml
```

### Resume a Failed Run

Nextflow caches completed tasks. To resume from the last successful step:

```bash
nextflow run main.nf -params-file params.yaml -resume
```

### Using Singularity Instead of Docker

```bash
nextflow run main.nf -params-file params.yaml -profile singularity
```

### Running on a Cluster

For SLURM:
```bash
nextflow run main.nf -params-file params.yaml -profile slurm
```

For SGE, LSF, or PBS, use `-profile sge`, `-profile lsf`, or `-profile pbs` respectively.

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--samplesheet` | Path to CSV samplesheet with columns: sampleID, fastq1, fastq2 |
| `--reference` | Path to reference genome FASTA file |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--fastp_qualified_quality` | `20` | Minimum quality for qualified bases |
| `--fastp_min_length` | `50` | Minimum read length after trimming |
| `--fastp_cut_mean_quality` | `20` | Mean quality threshold for sliding window |
| `--fastp_extra_args` | `""` | Additional fastp arguments |
| `--max_cpus` | `8` | Maximum CPUs per process |
| `--max_memory` | `32.GB` | Maximum memory per process |
| `--max_time` | `24.h` | Maximum time per process |

## Output Structure

```
results/
├── sample1/
│   ├── fastqc/
│   │   ├── raw/              # FastQC reports for raw reads
│   │   └── trimmed/          # FastQC reports for trimmed reads
│   ├── fastp/                # Trimming reports and trimmed FASTQ files
│   ├── alignment/            # BAM files (raw and sorted)
│   └── stats/                # Alignment statistics
├── sample2/
│   └── ...
├── multiqc/                  # Aggregated QC report
│   ├── multiqc_report.html
│   └── multiqc_data/
└── pipeline_info/            # Pipeline execution reports
    ├── timeline.html
    ├── report.html
    ├── trace.txt
    └── dag.svg
```

### Key Output Files

- **MultiQC Report**: `results/multiqc/multiqc_report.html` - Comprehensive QC overview
- **Sorted BAM**: `results/{sample}/alignment/{sample}.sorted.bam` - Final aligned reads
- **BAM Index**: `results/{sample}/alignment/{sample}.sorted.bam.bai` - BAM index file
- **Alignment Stats**: `results/{sample}/stats/` - flagstat, stats, and idxstats files
- **Pipeline Timeline**: `results/pipeline_info/timeline.html` - Execution timeline

## Customization

### Adjusting Trimming Parameters

For aggressive quality trimming:

```yaml
fastp_qualified_quality: 30
fastp_min_length: 75
fastp_cut_mean_quality: 30
fastp_extra_args: "--cut_front --cut_tail"
```

For lenient trimming (degraded samples):

```yaml
fastp_qualified_quality: 15
fastp_min_length: 35
fastp_cut_mean_quality: 15
```

### Adding Custom Adapters

```yaml
fastp_extra_args: "--adapter_fasta /path/to/adapters.fa"
```

### Enable Deduplication

```yaml
fastp_extra_args: "--dedup"
```

### Adjusting Resource Allocation

Edit `nextflow.config` to modify resources for specific processes:

```groovy
process {
    withName: BWA_MEM {
        cpus = 16
        memory = '32.GB'
        time = '12.h'
    }
}
```

## Troubleshooting

### Pipeline fails during reference indexing

If BWA index files don't exist, the pipeline will generate them. Ensure:
- Reference FASTA file path is correct
- You have write permissions in the reference directory

### Docker permission errors

If you encounter permission issues with Docker:

```bash
# Add your user to the docker group
sudo usermod -aG docker $USER
# Log out and back in for changes to take effect
```

### Out of memory errors

Reduce resource requirements in `params.yaml`:

```yaml
max_cpus: 4
max_memory: "16.GB"
```

Or edit specific process resources in `nextflow.config`.

### Sample validation errors

Ensure:
- Samplesheet is properly formatted (CSV with header: `sampleID,fastq1,fastq2`)
- All FASTQ file paths are absolute or relative to the pipeline directory
- All FASTQ files exist and are readable
- No spaces in sample IDs

## Quality Control Interpretation

### Key Metrics to Review in MultiQC

1. **Total Reads** - Check read counts are sufficient for your analysis
2. **Trimming Rate** - Percentage of reads surviving trimming (typically >80%)
3. **Mapping Rate** - Percentage of reads aligning (typically >80% for DNA-seq)
4. **Duplicate Rate** - PCR/optical duplicates (varies by library type)
5. **Insert Size** - Distribution should match expected library prep

### QC Thresholds (Recommended)

These are general guidelines; adjust based on your specific application:

- Minimum reads after trimming: 1M per sample
- Minimum mapping rate: 80%
- Maximum duplicate rate: 30%
- Per-base quality scores: >30 across most of the read

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316-319.
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890.
- **BWA**: Li, H. and Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 25, 1754-1760.
- **SAMtools**: Li, H., et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078-2079.
- **MultiQC**: Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32, 3047-3048.

## Support

For issues or questions:
1. Check the [Nextflow documentation](https://www.nextflow.io/docs/latest/)
2. Review error messages in `.nextflow.log`
3. Examine process-specific logs in `work/` directory

## License

This pipeline is released under the MIT License.

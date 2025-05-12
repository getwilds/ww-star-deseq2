# ww-star-deseq2
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This WILDS WDL workflow performs RNA-seq analysis using STAR's two-pass alignment methodology and DESeq2 differential expression analysis. It is intended to be a straightforward demonstration of an RNA sequencing pipeline within the context of the WILDS ecosystem.

## Overview

The workflow performs the following key steps:
1. Optional automatic reference genome download (if not provided)
2. STAR index building
3. STAR two-pass alignment for each sample
4. RNA-SeQC quality control analysis
5. Combining gene count matrices
6. DESeq2 differential expression analysis with visualization

## Features

- Two-pass STAR alignment for improved splice junction detection
- Automatic reference genome download (optional)
- Quality control with RNA-SeQC
- Differential expression analysis with DESeq2
- Visualization of results (PCA, volcano plots, heatmaps)
- Compatible with the WILDS workflow ecosystem

## Usage

### Requirements

- [Cromwell](https://cromwell.readthedocs.io/), [MiniWDL](https://github.com/chanzuckerberg/miniwdl), or another WDL-compatible workflow executor
- Docker/Apptainer (the workflow uses WILDS Docker containers)

### Basic Usage

1. Create an inputs JSON file with your sample information:

```json
{
  "star_deseq2.samples": [
    {
      "name": "sample1",
      "r1": "/path/to/sample1_1.fastq.gz",
      "r2": "/path/to/sample1_2.fastq.gz",
      "condition": "treatment"
    },
    {
      "name": "sample2",
      "r1": "/path/to/sample2_1.fastq.gz",
      "r2": "/path/to/sample2_2.fastq.gz",
      "condition": "control"
    }
  ],
  "star_deseq2.reference_level": "control",
  "star_deseq2.contrast": "condition,treatment,control"
}
```

2. Run the workflow using your preferred WDL executor:

```bash
# Cromwell
java -jar cromwell.jar run ww-star-deseq2.wdl --inputs ww-star-deseq2-inputs.json --options ww-star-deseq2-options.json

# miniWDL
miniwdl run ww-star-deseq2.wdl -i ww-star-deseq2-inputs.json
```

### Integration with SRA Downloads

The workflow pairs well with the [ww-sra](https://github.com/getwilds/ww-sra) workflow for downloading data from NCBI's Sequence Read Archive.

### Detailed Options

The workflow accepts the following inputs:

| Parameter | Description | Type | Required? | Default |
|-----------|-------------|------|-----------|---------|
| `samples` | Array of sample information objects | Array[SampleInfo] | Yes | - |
| `reference_genome` | Reference genome information | RefGenome | No | GRCh38.p14 (auto-downloaded) |
| `reference_level` | Reference level for DESeq2 | String | No | "" |
| `contrast` | Contrast string for DESeq2 | String | No | "" |

#### SampleInfo Structure

Each entry in the `samples` array should contain:
- `name`: Sample identifier
- `r1`: Path to R1 FASTQ file
- `r2`: Path to R2 FASTQ file
- `condition`: Group/condition for differential expression

#### RefGenome Structure (optional)

If provided, the `reference_genome` should contain:
- `name`: Reference genome name
- `fasta`: Path to reference FASTA file
- `gtf`: Path to reference GTF file

### Output Files

The workflow produces the following outputs:

| Output | Description |
|--------|-------------|
| `star_bam` | Aligned BAM files for each sample |
| `star_bai` | BAM indexes for each sample |
| `star_gene_counts` | Raw gene counts for each sample |
| `star_log_final` | STAR final logs |
| `star_log_progress` | STAR progress logs |
| `star_log` | STAR main logs |
| `star_sj` | STAR splice junction files |
| `rnaseqc_metrics` | RNA-SeQC quality metrics |
| `combined_counts_matrix` | Combined gene counts matrix |
| `sample_metadata` | Sample metadata table |
| `deseq2_all_results` | Complete DESeq2 results |
| `deseq2_significant_results` | Filtered significant DESeq2 results |
| `deseq2_normalized_counts` | DESeq2 normalized counts |
| `deseq2_pca_plot` | PCA plot of samples |
| `deseq2_volcano_plot` | Volcano plot of differential expression |
| `deseq2_heatmap` | Heatmap of differentially expressed genes |

## For Fred Hutch Users

For Fred Hutch users, we recommend using [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster. To do this:

1. Start by either cloning or downloading a copy of this repository to your local machine.
    - Cloning: `git clone https://github.com/getwilds/ww-star-deseq2.git`
    - Downloading: Click the green "Code" button in the top right corner, then click "Download ZIP".
2. Update `ww-star-deseq2-inputs.json` with your sample names, FASTQ file paths, and conditions.
3. Update `ww-star-deseq2-options.json` with your preferred location for output data to be saved to (`final_workflow_outputs_dir`).
4. Submit the WDL file along with your custom json's to the Fred Hutch cluster via PROOF by following our [SciWiki documentation](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/).

Additional Notes:
- Keep in mind that all file paths in the jsons must be visible to the Fred Hutch cluster, e.g. `/fh/fast/`, AWS S3 bucket. Input file paths on your local machine won't work in PROOF.
- Specific reference genome files can be provided as inputs, but if none are provided, the workflow will automatically download a GRCh38 reference genome and use that. For the first go-around, we recommend starting with the default reference files.
- To avoid duplication of reference genome data, we highly recommend executing this workflow with call caching enabled in the options json (`write_to_cache`, `read_from_cache`, already set to `true` here).

## Docker Containers

This workflow uses the following Docker containers from the [WILDS Docker Library](https://github.com/getwilds/wilds-docker-library):

- `getwilds/star:2.7.6a` - For STAR alignment
- `getwilds/rnaseqc:2.4.2` - For RNA-SeQC quality control
- `getwilds/gtf-smash:v8` - For reference genome download
- `getwilds/combine-counts:0.1.0` - For combining count matrices
- `getwilds/deseq2:1.40.2` - For differential expression analysis

All containers are available on both DockerHub and GitHub Container Registry.

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on our [issue tracker](https://github.com/getwilds/ww-star-deseq2/issues).

## Contributing

If you would like to contribute to this WILDS WDL workflow, see our [contribution guidelines](.github/CONTRIBUTING.md) as well as our [WILDS Contributor Guide](https://getwilds.org/guide/) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

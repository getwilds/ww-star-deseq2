
# ww-star-deseq2
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)

This WILDS WDL workflow performs alignment using the [two-pass methodology](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) of [STAR](https://github.com/alexdobin/STAR) and subsequently analyzes that alignment via [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). It is intended to be a relatively straightforward demonstration of an RNA sequencing pipeline within the context of the WILDS ecosystem.

## Basic Usage

For Fred Hutch users that are new to WDL, we recommend using [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster, as it simplifies interaction with Cromwell and provides a user-friendly front-end for job submission and tracking. To do this:

1. Start by cloning/downloading a copy of this repository to your local machine (click "Code", then click "Download ZIP").
2. Update [`ww-star-deseq2-inputs.json`](https://github.com/getwilds/ww-star-deseq2/blob/main/ww-star-deseq2-inputs.json) with your sample names (`omics_sample_name`) and FASTQ file paths (`R1` and `R2`).
3. Update [`ww-star-deseq2-options.json`](https://github.com/getwilds/ww-star-deseq2/blob/main/ww-star-deseq2-inputs.json) with your preferred location for output data to be saved to (see `final_workflow_outputs_dir` field).
4. Submit the WDL file along with your custom json's to the Fred Hutch cluster via PROOF by following our [SciWiki documentation](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/).

Additional Notes:
- Keep in mind that all file paths in the jsons must be visible to the Fred Hutch cluster, e.g. `/fh/fast/`, AWS S3 bucket. Input file paths on your local machine won't work in PROOF.
- Specific reference genome files can be provided as inputs, but if none are provided, the workflow will automatically download a GRCh38 reference genome and use that. For the first go-around, we recommend starting with the default reference files.

## Advanced Usage

For users outside of Fred Hutch or more advanced users who would like to run the workflow locally, command line execution is relatively straightforward: 
```
java -jar cromwell-86.jar run ww-star-deseq2.wdl --inputs ww-star-deseq2-inputs.json --options ww-star-deseq2-options.json
```
Although Cromwell is demonstrated here, this pipeline is not specific to Cromwell and can be run using whichever WDL execution method you prefer ([miniwdl](https://github.com/chanzuckerberg/miniwdl), [Terra](https://terra.bio/), [HealthOmics](https://docs.aws.amazon.com/omics/latest/dev/workflows.html), etc.).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on our [issue tracker](https://github.com/getwilds/ww-star-deseq2/issues).

## Contributing

If you would like to contribute to this WILDS WDL workflow, see our [contribution guidelines](.github/CONTRIBUTING.md) as well out our [WILDS Contributor Guide](https://getwilds.org/guide/) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.

# autocycler-nf

A pipeline for running [rrwick/Autocycler](https://github.com/rrwick/Autocycler). Based on Ryan Wick's [autocycler_full.sh](https://github.com/rrwick/Autocycler/tree/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick) script.

## Development Status

In active development. Not ready for use.

## Usage

### Basic usage

```
nextflow run BCCDC-PHL/autocycler-nf \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input_long /path/to/your/fastqs \
  --outdir /path/to/your/outputs
```

By default, the pipeline will divide the input fastqs into 4 subsamples and assemble each subsample using the following assemblers:

```
raven
myloasm
miniasm
flye
metamdbg
necat
nextdenovo
plassembler
canu
```

...resulting in 4 x 9 = 36 assemblies for each sample, which are combined into a final 'consensus assembly'.


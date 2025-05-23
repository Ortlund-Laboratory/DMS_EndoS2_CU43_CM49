# Example of Workflow to Calculate Binding Scores

Fc variant binding to CM49 was investigated, and reference, tighter-binding and weaker-binding populations collected for sequencing. These sequences can be processed to give scores for both tighter and weaker binders.

## Input Files Required

[SnakeFile](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/Snakefile)<br>
Gives overall instructions for the `snakemake` workflow.<br>
[config.yaml](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Markdown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/data/barcode_runs.csv)<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/data/processed_ccs.csv)<br>
Processed PacBio CCSs, generated from our [PacBio_Library_Sequencing](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main/PacBio_Library_Sequencing) routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/data/wildtype_sequence.fasta)<br>
Fc wildtype sequence.<br>

### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/data/barcode_runs.csv). These files are too large to be contained in GitHub, and so are found, respectively, at:

**NCBI**<br>
Overall BioProject: [PRJNA1263835](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1263835)

**EndoS2**<br>
BioSample: [SAMN48540912](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN48540912)

**CM49**<br>
BioSample: [SAMN48540122](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN48540122)

**CU43**<br>
BioSample: [SAMN48540644](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN48540644)

**Reference for EndoS2 experiments**<br>
BioSample: [SAMN48541669](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN48541669)

**Reference for CM49 and CU43 experiments**<br>
BioSample: [SAMN48541626](https://www.ncbi.nlm.nih.gov/biosample/?term=SAMN48541626)

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Key Output

[results/counts/barcode_fates.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/results/counts/barcode_fates.csv)<br>
Tally of barcodes classified and filtered according to quality.<br>
[results/counts/variant_counts.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/results/counts/variant_counts.csv)<br>
Tally of individual barcode counts for each sample.<br>

## Binding Score Generation & Data Visualization

Go to [scores_and_visualization](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization) for tight-binding and/or weaker-binding score calculation and heatmap generation.

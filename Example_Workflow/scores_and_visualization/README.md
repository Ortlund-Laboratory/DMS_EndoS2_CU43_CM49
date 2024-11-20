# File Prep

[variant_counts.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/results/counts/variant_counts.csv) from results/counts needs to be manipulated.

In our example, by referring to [barcode_runs.csv](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/blob/main/example_enrichment_escape/data/barcode_runs.csv), we can see that three experiments pertinent to CM49 were run: a reference sample (exp01-none-0-reference), a sample for tighter binding with CM49 (exp04-CM49top-2000-escape) and a sample for weaker binding with CM49 (exp05-CM49bottom-2000-escape). 

**NOTE**: Though the above samples have escape in the title, this was due to naming conventions in some of the Jupyter scripts, which requires samples to be classified as either 'reference' or 'escape'. In our naming system, samples which have 'top' in their titles denote tighter binding populations, and samples with 'bottom' in their titles denote weaker binding populations.

For compatibility with our R scripts, we require the information in [variant_counts.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/results/counts/variant_counts.csv) to be split into three-column (barcode,mutation,count) files. This is very simple and you can write a script to do this if you like.

First, remove all variants where the number of amino acid mutations is 0 (i.e. is still wildtype) or greater than 1 (multiple point mutations).

```
awk -F',' '!($10!=1)' variant_counts.csv > tmp.csv && mv tmp.csv variant_counts.csv
```
Next step is to separate into reference, tighter-binding (enrichment) and weaker-binding (escape) files.
```
sed -n '/exp01-none-0-reference/p' variant_counts.csv > ref_variant_counts.txt
sed -n '/exp04-CM49top-2000-escape/p' variant_counts.csv > enrich_variant_counts.txt
sed -n '/exp05-CM49bottom-2000-escape/p' variant_counts.csv > escape_variant_counts.txt
```
Then reformat each file so that they just contain tab-separated columns for barcode, mutation and count.
```
awk '{print $4, $8, $5}' FS="," OFS="\t" ref_variant_counts.txt > tmp.txt && mv tmp.txt ref_variant_counts.txt
awk '{print $4, $8, $5}' FS="," OFS="\t" enrich_variant_counts.txt > tmp.txt && mv tmp.txt enrich_variant_counts.txt
awk '{print $4, $8, $5}' FS="," OFS="\t" escape_variant_counts.txt > tmp.txt && mv tmp.txt escape_variant_counts.txt
```
Finally, add tab-separated headers (barcode  mutation  count) for each of these files.

Now that we have the input files formatted correctly, we are ready to calculate binding scores and to visualize the data. We run these steps in separate tighter and weaker subdirectories.


## Input Files Required

[BarcodeMapping_CM49bottom.R](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/weaker/BarcodeMapping_CM49bottom.R)<br>
R script to calculate scores from reference and escape counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/weaker/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[escape_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/weaker/escape_variant_counts.txt)<br>
Barcode counts for the escape sample.<br>
[Fc_prot.fasta](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/weaker/Fc_prot.fasta)<br>
Amino acid sequence for WT Fc. This is required to complete the heatmap.

## Workflow

```
rstudio BarcodeMapping_CM49bottom.R
```

## Key Output

[CM49bottom_Weaker_fractions.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/weaker/output/CM49bottom_Weaker_fractions.csv)<br>
Log of each mutation and its associated weaker-binding score.<br>
[CM49bottom_WeakerFraction_heatmap01.png](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/weaker/output/CM49bottom_WeakerFraction_heatmap01.png)<br>
Heatmap of weaker-binding scores for all single-point Fc variants.<br>

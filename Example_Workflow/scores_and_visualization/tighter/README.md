## Input Files Required

[BarcodeMapping_CM49top.R](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/tighter/BarcodeMapping_CM49top.R)<br>
R script to calculate scores from reference and enrichment counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/tighter/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[enrich_variant_counts.txt](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/tighter/enrich_variant_counts.txt)<br>
Barcode counts for the enrichment sample.<br>
[Fc_prot.fasta](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/tighter/Fc_prot.fasta)<br>
Amino acid sequence for WT Fc. This is required to complete the heatmap.

## Workflow

```
rstudio BarcodeMapping_CM49top.R
```

## Key Output

[CM49top_Tighter_fractions.csv](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/tighter/output/CM49top_Tighter_fractions.csv)<br>
Log of each mutation and its associated tighter-binding score.<br>
[CM49top_TighterFraction_heatmap01.png](https://github.com/Ortlund-Laboratory/DMS_EndoS2_CU43_CM49/blob/main/Example_Workflow/scores_and_visualization/tighter/output/CM49top_TighterFraction_heatmap01.png)<br>
Heatmap of tighter-binding scores for all single-point Fc variants.<br>

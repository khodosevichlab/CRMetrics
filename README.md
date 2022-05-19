# CRMetrics
Summarization of Cell Ranger count metrics. 

As example data, you can use /data/PD-MSA_lentiform_nucleus/counts_premrna

The metrics are stored in <SAMPLE>/outs/metrics_summary.csv

Notice one sample is marked "FAIL" and should not be included. Consider including an exclude parameter.

As Laura mentioned, it would be nice to be able to show distribution plots of e.g. total UMIs per cell/sample. Therefore, consider including an option to load count matrices (slower) than to just load metrics summaries.

Please check inspiration.Rmd.
